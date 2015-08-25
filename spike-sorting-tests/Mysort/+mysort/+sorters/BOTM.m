classdef BOTM < mysort.sorters.OnlineSpikeSorterInterface & ...
                mysort.sorters.FilterBasedSpikeSorterInterface
    
    properties (SetAccess=private)
        init_counter
    end
    properties
        E
        D
        lastChunkSortingD
        channelSets
        MAHA
        CONF_up
        Tf_up
        blockLen
            
        noisePrior
        threshold
        priors
        Dshifts
        templateAlignmentShifts
         
        spikeEpochs
        
        SICvaluesForDebugging
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = BOTM(Covest, T, varargin)
            if max(size(T)) == 1
                % be downward compatible with old calling syntax where T is
                % in vector concatenated form and Tf was supplied first
                Tf = T;
                T = varargin{1};
                varargin(1) = [];
            else % T is in tensor form, Tf x nC x nWfs
                Tf = size(T,1);
                T = mysort.wf.t2v(T);
            end            
            
            self = self@mysort.sorters.OnlineSpikeSorterInterface(varargin{:});
            self = self@mysort.sorters.FilterBasedSpikeSorterInterface(Covest, Tf, T, varargin{:});
            %self = self@mysort.spiketrain.SpikeSortingContainer('BOTM', []);
            self.P.spikePrior = 0.00001;
            self.P.upsample = 3;
            self.P.adapt = false;
            self.P.adaptOnInit = false;
            self.P.templMaxIdx = [];
            self.P.max_num_of_channels_per_template = 10;
            self.P.norm_constraint = false;
            self.P.verbose = self.LEVEL_PROCESSING;
            self.P.useSIC = true;
            self.P.storeSICvaluesForDebugging = 0;
            self.P.minEpochLength = round(.5*Tf);
            self.P = mysort.util.parseInputs(self.P, varargin, 'error');
            
            if ~self.P.useSIC
                % There is no need to upsample if we do not use SIC
                self.P.upsample = 1;
            end
            
            if isempty(self.P.templMaxIdx)
                TT = mysort.wf.v2m(T, size(T,2)/Tf);
                [i j] = mysort.util.matrixArgMax(abs(TT));
                self.P.templMaxIdx = j;
            end
            
            self.E = [];  
            self.D = [];
            self.lastChunkSortingD = [];
            self.MAHA = [];   
            self.channelSets = {};
            self.CONF_up = [];  
            self.Dshifts = [];
            self.Tf_up = [];
            self.blockLen = [];
            
            self.noisePrior = [];
            self.threshold = [];
            self.priors = [];
            self.activeSpikeSortingIdx = 1;
            if self.P.adaptOnInit
                self.adapt(1);
            end
        end
        
        %%% ------------------------------------------------------
        function [sorting filterOut] = sortMatrix(self, CHUNK)
            filterOut = self.getFilterOutputs_(CHUNK);
            self.debugout('Postprocessing...', self.LEVEL_PROCESSING);
            sorting = self.Postprocessing();
        end
        %%% ------------------------------------------------------
        function b = hasSpikeSorting(self)
            b = true;
        end
        %%% ------------------------------------------------------
        function S = getActiveSpikeSorting(self)
            S = self;
        end
        %%% ------------------------------------------------------
        function T = getTemplateWaveforms(self)
            T = mysort.wf.v2t(self.T, size(self.DH,2));
        end
        %------------------------------------------------------------------
        function cutleft = getTemplateCutLeft(self)
            cutleft = ceil(self.Tf/2)-1;
        end
        %------------------------------------------------------------------
        function sorting = getGdf(self, t1, t2, unitNames)
            assert(nargin == 3, 'This function is only implemented for 3 parameters!');
             [sorting filterOut] = sortMatrix(self, self.DH(t1:t2,:)');
             if ~isempty(sorting)
                sorting(:,2) = sorting(:,2) + t1-1;
             end
        end        
        
        
        %%% ------------------------------------------------------
        function fout = getFilterOutputs_(self, X, filterindex)
            self.debugout('Applying filter...', self.LEVEL_PROCESSING);
            self.Y = mysort.util.applyMcFilters(X, self.F);

            self.debugout('Calculating Discriminant Function...', self.LEVEL_PROCESSING);
            [self.D self.Dshifts] = mysort.util.calculateDiscriminantFunctions(self.Y, self.CONF, self.priors);  
            
            if 0
                ax = [];
                mysort.plot.figure('w', 800, 'h', 800);
                ax(1) = subplot(4,1,1)
                hold on
                for i=1:size(self.D,1)
                    plot(self.D(i,:), 'color', mysort.plot.vectorColor(i), 'linewidth', 2);
                end
                title('BOTM discrm func')
                
                [Da alphas] = mysort.util.calculateDiscriminantFunctionsWithAlpha(self.Y, self.CONF, self.priors);  
                Palpha = evpdf(alphas, 1, .2);
                Palpha(alphas<.1) = 0;
                logPalpha = log(Palpha);
                logPalpha(logPalpha < -2) = nan;
                ax(2) = subplot(5,1,2);
                hold on
                for i=1:size(self.D,1)
                    plot(alphas(i,:), 'color', mysort.plot.vectorColor(i), 'linewidth', 2);
                end
                title('optimal alphas (scaling of templates)')
                set(ax(2), 'ylim', [0 1.5]);
                
                ax(3) = subplot(5,1,3);
                hold on
                for i=1:size(self.D,1)
                    plot(Palpha(i,:), 'color', mysort.plot.vectorColor(i), 'linewidth', 2);
                end
                title('p(alpha)')
                set(ax(3), 'ylim', [0 1.5]);
                
                ax(4) = subplot(5,1,4)
                hold on
                for i=1:size(self.D,1)
                    plot(Da(i,:), 'color', mysort.plot.vectorColor(i), 'linewidth', 2);
                end
                title('Discr wo alpha prior')
                
                ax(5) = subplot(5,1,5)
                hold on
                for i=1:size(self.D,1)
                    plot(Da(i,:)+logPalpha(i,:), 'color', mysort.plot.vectorColor(i), 'linewidth', 2);
                end
                linkaxes(ax, 'x');
                [maxis maxChan] = max(Da+logPalpha,[], 1);
                [pks locs] = findpeaks(maxis,'MINPEAKDISTANCE', 2);
                ids = maxChan(locs);
                for i=1:size(self.D,1)
                    plot(locs(ids==i), pks(ids==i), 'o', 'color', mysort.plot.vectorColor(i), 'linewidth', 2);
                end
                title('Discr funcs')
            end
            
            if nargout > 0
                fout = self.D;
                if nargin == 3
                    fout = fout(filterindex, :);
                end
            end     
            if size(X, 2) < size(fout, 2)
                fout = fout(:, end-size(X, 2)+1:end);
            end
        end
        
        %%% ------------------------------------------------------
        %%% -----------------PROCESSING FUNCTIONS-----------------
        %%% ------------------------------------------------------  
        
        %%% ------------------------------------------------------        
        function adapt(self, initFlag)
            % In this version adapt only the first time
            if nargin == 1
                if (self.P.adaptOnInit || self.currentChunk > 1) && ~self.P.adapt
                    return
                end
            end
            % "rescue" the original templates
            self.templates = self.T;
            nC = size(self.Covest.CCol,2);
            
            if self.currentChunk > 1
                self.debugout('Adapting Templates...', self.LEVEL_PROCESSING); 
                maxPastSample = max(1, self.chunk_end - 20*self.DH.getSamplesPerSecond());
                classes = self.getSpikeClasses('start', maxPastSample ,'stopp', self.chunk_end);
                epochs = self.getSpikeEpochs('start', maxPastSample ,'stopp', self.chunk_end);
                X = self.DH(maxPastSample:self.chunk_end, :)';
                for t = 1:size(self.T,1)
                    myNonOverlappingSpikes = mysort.epoch.removeOverlapping(epochs(classes==t,:), epochs(classes~=t,:));
                    weight = 1;
                    if size(myNonOverlappingSpikes,1) > 250
                        myNonOverlappingSpikes = myNonOverlappingSpikes(end-250+1:end,:);
%                     elseif size(myNonOverlappingSpikes,1)>0
%                         weight = min(1, size(myNonOverlappingSpikes,1)/10);
                    elseif isempty(myNonOverlappingSpikes)
                        continue
                    end
                    spikes = mysort.epoch.extractWaveform(X, myNonOverlappingSpikes - maxPastSample +1);
                    UP = 3;
                    spikesUP = mysort.util.resampleTensor(mysort.util.m2t(spikes, nC), UP,1);
                    spikesUP = mysort.wf.t2m(spikesUP);
                    [tau spikesUP] = mysort.util.alignWaveformsOnMaxOrMin(spikesUP, nC, 'maxIdx', UP*self.P.templMaxIdx);
                    t_up = mean(spikesUP, 1);
                    self.T(t,:) = (1-weight) * self.T(t,:) + ...
                                     weight  * mysort.wf.m2v(mysort.util.resampleMC(mysort.util.v2m(t_up, nC), 1, UP));
                end
            end
%             [self.templateAlignmentShifts self.T] = mysort.util.alignWaveformsOnMaxOrMin(self.T, nC, 'maxIdx', self.P.templMaxIdx);
            tT = mysort.wf.v2t(self.T, nC);
            [ self.T self.templateAlignmentShifts] = mysort.wf.tAlignOnCorrelation(tT, 'trunc', 1, 'absMax', 1);
            self.T = mysort.wf.t2v(self.T);
%             self.T = -self.T;
            % set all channels to exactly zero, if there is not enough
            % energy on the whole channel in that template
            % maxAmpl = max(abs(T(:)));
            % Calculate optimally matched Filters
            self.debugout('Calculating Filters...', self.LEVEL_PROCESSING);            
            amplThresh = .1;
            self.channelSets = {};
            for t=1:size(self.T,1)
                % get one template
                templ = mysort.wf.v2m(self.T(t,:), nC);
                
                % compute the mahalanobis energies on every channel
                % individually
                mahalanobis_energy = zeros(nC, 2);
                for chan = 1:size(templ,1)
                    nC_red = 1;
                    t_red = templ(chan,:);
                    if ~any(t_red)
                        mahalanobis_energy(chan, 2) = 0;
                    else
                        t_red_ = mysort.util.embedTime2embedChan(t_red, nC_red);
                        CCol_red = mysort.noise.ccolSubChanIdx(self.Covest.CCol, chan, self.Tf-1);
                        CCol_red(1:nC_red, 1:nC_red) = CCol_red(1:nC_red, 1:nC_red) + 0*diag(ones(1,nC_red));
                        f_red_ = matlabfilecentral.block_levinson(t_red_(:), CCol_red);
                        f_red = mysort.util.embedChan2embedTime(f_red_(:), nC_red);                    
                        mahalanobis_energy(chan, 2) = sqrt(t_red*f_red/self.Tf);
                    end
                end
                mahalanobis_energy(:,1) = [1:nC]';
                mahalanobis_energy = sortrows(mahalanobis_energy, -2);

                % set the template to zero on all channels with too low SNR
                
                % take highes channel for sure:
                self.channelSets{t} = mahalanobis_energy(1,1);
                templ_reduced = zeros(size(templ));
                templ_reduced(mahalanobis_energy(1,1), :) = templ(mahalanobis_energy(1,1), :);
                energies = mahalanobis_energy(1,2);
                
                % take up to self.P.max_num_of_channels_per_template-1 more channels
                for chan = 2:min(self.P.max_num_of_channels_per_template, size(mahalanobis_energy,1))
                    chanIdx = mahalanobis_energy(chan,1);
                    if mahalanobis_energy(chan,2) > amplThresh  || ...
                       mahalanobis_energy(chan,2) > sum(mahalanobis_energy(:,2))/20
                        templ_reduced(chanIdx,:) = templ(chanIdx,:);
                        self.channelSets{t} = [self.channelSets{t} chanIdx];
                        energies = [energies mahalanobis_energy(chan,2)];
                    else
                        break
                    end
                end
                fprintf('Using %d channels for template %d (energies: ', length(self.channelSets{t}), t);
                fprintf('%.3f ', energies);
                fprintf(')\n');

%                 clf
%                 plot(mahalanobis_energy(:,1), mahalanobis_energy(:,2), 'kx')
%                 hold on
%                 plot([0 size(templ,1)], [amplThresh amplThresh], 'g:');
%                 plot([0 size(templ,1)], sum(mahalanobis_energy(:,2))/20 *[1 1], 'r:');
%                 plot(self.channelSets{t}, energies, 'xm', 'markersize', 15, 'linewidth', 3);
%                 
                % build the "reduced" template:
                nC_red = length(self.channelSets{t});
                t_red = mysort.wf.m2v(templ(self.channelSets{t},:));
                if ~any(t_red)
                    f_red = mysort.wf.v2m(t_red, nC_red);
                else
                    t_red_ = mysort.util.embedTime2embedChan(t_red, nC_red);
                    CCol_red = mysort.noise.ccolSubChanIdx(self.Covest.CCol, self.channelSets{t}, self.Tf-1);
                    CCol_red(1:nC_red, 1:nC_red) = CCol_red(1:nC_red, 1:nC_red) + 0*diag(ones(1,nC_red));
                    f_red_ = matlabfilecentral.block_levinson(t_red_(:), CCol_red)';
                    f_red = mysort.util.embedChan2embedTime(f_red_, nC_red);
                    t_red = mysort.wf.v2m(t_red, nC_red);
                    f_red = mysort.wf.v2m(f_red, nC_red);
                    if 0
                        figure;
                        subplot(4,1,1)
                        plot(t_red');
                        subplot(4,1,2)
                        plot(f_red');
                        subplot(4,1,3)
                        xc = mysort.util.mcfilt(t_red, f_red);
                        plot(xc);
                        subplot(4,1,4)
                        xc = mysort.util.mcfilt(t_red+40*randn(size(t_red)), f_red);
                        plot(xc);                        
                    end                
                end                    
%                 % get "reduced" covariance matrix
%                 Cinv_red = self.NE.inv(self.Tf, self.channelSets{t});
%                 % compute "reduced" filter
%                 F_red = T_red*Cinv_red;
%                 f_red = mysort.wf.v2m(F_red, length(self.channelSets{t}));
                
                % build the full filter
                f = zeros(nC, self.Tf);
                for chan =1:size(templ,1)
                    reduced_chan = find(self.channelSets{t}==chan);
                    if ~isempty(reduced_chan)
                        f(chan,:) = f_red(reduced_chan,:);
                    end
                end
                
                self.T(t,:) = mysort.wf.m2v(templ);
                self.F(t,:) = mysort.wf.m2v(f);
            end

            self.debugout('Done.', self.LEVEL_PROCESSING);            
            self.nF = size(self.F,1);
            % Calculate Priors            
            self.noisePrior = 1-self.nF * self.P.spikePrior; 
            self.priors = repmat(self.P.spikePrior, self.nF, 1); 
            self.threshold = log(self.noisePrior);
            % Calculate confusion matrix
            self.CONF    = mysort.util.calculateXIvsF(self.T, self.F, nC, 0);
            self.CONF_up = mysort.wf.tResample(self.CONF, self.P.upsample, 1);
%             self.CONF_up = mysort.util.resampleTensor(self.CONF, self.P.upsample, 1);

            self.Tf_up = self.Tf * self.P.upsample;
            self.blockLen = ceil(self.Tf_up/4);  
            %mysort.plot.XIvsF(mysort.util.m2t(self.T, 4), mysort.util.m2t(self.F, 4),'XIvsF', self.CONF);
        end

        %%% ------------------------------------------------------   
        function chunk_sorting = Postprocessing(self)
            % Does the postprocessing of the filter calculation for this
            % CHUNK. D are the Discriminant functions of the current chunk,
            % CHUNK is the original data with overlapping regions to the
            % left and right CHUNK.  
            [mPr c] = max(self.D,[],1);
            spikeEpochs = mysort.epoch.fromBinaryVectorMinLen(...
                                    mPr>self.threshold, self.P.minEpochLength);
            chunk_sorting = self.processSpikeEpochs(spikeEpochs);
            if ~isempty(chunk_sorting)
                chunk_sorting = sortrows(chunk_sorting,2);
            end
        end
        
        %%% ------------------------------------------------------   
        function chunk_gdf = processSpikeEpochs(self, epochs)
            chunk_gdf = [];
            if isempty(epochs); return; end
            nE = size(epochs,1);      
            %  figure; plot(self.D')
            if self.P.storeSICvaluesForDebugging
                self.SICvaluesForDebugging = [];
            end
            for e=1:nE
                if self.P.upsample>1
                    Yup = [];
                    % UPSAMPLE THE FILTER OUTPUTS NOT THE DISCRIMINANT
                    % FUNCTIONS !!! D will be not zero mean and you will get
                    % resampling artifacts at the right end!
                    
                    % Cut a bit wider at the edges and throw the edge
                    % resampling artefacts away
                    o1 = epochs(e,1)-max(1, epochs(e,1)-10);
                    o2 = min(size(self.Y,2), epochs(e,2)+10) - epochs(e,2);
                    e1 = epochs(e,1)-o1;
                    e2 = epochs(e,2)+o2;
                    Yup(1:self.nF,:) = mysort.util.resampleMC(...
                        self.Y(1:self.nF,e1:e2), self.P.upsample, 1);
                    Yup = Yup(:, o1*self.P.upsample+1:end-o2*self.P.upsample);
                else
                    e1 = epochs(e,1);
                    e2 = epochs(e,2);
                    Yup = self.Y(:,e1:e2);
                end
                [epoch_gdf, SICvals] = self.processEpoch(Yup, self.CONF_up, self.Tf_up, self.priors, epochs(e,:));
                if ~isempty(epoch_gdf)
                    epochOffset = epochs(e,1)-ceil(self.Tf/2)-1;
                    epoch_gdf(:,2) = round(epoch_gdf(:,2)./self.P.upsample) + epochOffset;
                    chunk_gdf = [chunk_gdf; epoch_gdf];
                    SICvals.epoch = [e1 e2];
%                     if self.P.storeSICvaluesForDebugging && e==1
%                         self.SICvaluesForDebugging = SICvals;
%                     else
%                         self.SICvaluesForDebugging(e) = SICvals;
%                     end
                end
            end
        end

        %%% ------------------------------------------------------   
        function [gdf SICvals] = processEpoch(self, Y, M, Tf, p, epoch)
            % Detect and remove Spikes as long as discrimant functions
            % are larger than the noise prior
            % - uses SIC
            D = mysort.util.calculateDiscriminantFunctions(Y, M, p);
            [maxD maxClasses] = max(D,[],1); 
            gdf = [];
            count = 0;
            debug = 0;
            if debug
                close all
                %mysort.plot.XIvsF(randn(3,Tf-2), randn(3,Tf-2), 'XIvsF', self.CONF_up)
            end
            
            SICvals = struct();

            if ~self.P.useSIC
                [maxd t] = max(maxD,[],2);
                maxC = maxClasses(t);
                offset = t - (Tf - (self.P.upsample-1));
                % correct for any shift that was done to the waveform
                % in the initialization
                tau = self.P.upsample*self.templateAlignmentShifts(maxC);
                gdf = [gdf; [maxC t+self.P.upsample+1+tau]];               
                if isempty(gdf)
                    fprintf('No spikes in Epoch\n');
                end
                return
            end

            subtractor = [];
            D_old = [];
            while any(maxD>self.threshold)
                count = count +1;
%                 if count >= 2
%                     debug = 1;                    
%                 end
                if count > 10
                    %error('too many!');
%                     warning('detected more than 10 spikes detected in this epoch!');
                end
                [maxd t] = max(maxD,[],2);
                maxC = maxClasses(t);
                offset = t - (Tf - (self.P.upsample-1));
                subtractor_old = subtractor;
                subtractor = zeros(size(D));
                for f=1:self.nF
                    sub = squeeze(M(:,f,maxC))';
                    if f == maxC
                        % Block double self detections?
                        subtractor(f,:) = -inf;
                    else                    
                        subtractor(f,:) = mysort.util.shiftSubtract(zeros(1,size(D,2)),...
                                        sub, offset, false);
                    end
        
                end
                
                % Check if the subtraction of the expected response is
                % indeed decreasing the norm of the filter outputs.
                % Does not work for MEA data yes. Set "norm_constraint" to
                % "false"
                if ~self.P.norm_constraint || ...
                        norm(Y) > norm(Y+subtractor)
                    if debug
                        fig = mysort.plot.figure();
                        ah = mysort.plot.subplot([3,1]);
                        axes(ah(1));
                        title('X');
                        Margin = 200*self.P.upsample;
                        Margin_up = 200*self.P.upsample;
                        Tf_up = Tf*self.P.upsample;
                        
%                         xrange = Margin_up:self.P.upsample:epoch(2)-epoch(1)+Margin_up;
                        x = self.DH(epoch(1)-Margin_up:epoch(2)+Margin_up,:)';
                        plot(x,...
                            'k', 'linewidth', 3);
                        
                        axes(ah(2));
                        title(['D' num2str(count)]);
                        plot([1 size(D,2)], [0 0], ':k');
                        hold on                                        
                        if ~isempty(subtractor_old)
                            plot(D_old','--', 'linewidth',3);
                        end   
                        plot(D','linewidth',3);
%                         spacer = mysort.plot.mc(D,'plotZeroLine',1,...
%                             'linewidth',3,'figure',0,'color',{'k'}, 'figure', 0);
        
                    end
                    D_old = D;
                    D = D + subtractor + log(p(maxC));
                    if debug
                        %mysort.plot.mc(D, 'spacer', spacer,'linewidth',2 , 'color', {'r'},'figure',0);
                        plot(D',':', 'linewidth',2);
                        axes(ah(3));                        
                        title('Subtractor');
                        %mysort.plot.mc(subtractor + log(self.priors(maxC)), 'spacer', spacer,'linewidth',2,'figure',0,'color', {'g'});
                        plot((subtractor + log(self.priors(maxC)))', 'linewidth',2);
                        linkaxes(ah, 'x');
                        axis tight
                    end
                    % correct for any shift that was done to the waveform
                    % in the initialization
                    tau = self.P.upsample*self.templateAlignmentShifts(maxC);
                    gdf = [gdf; [maxC t+self.P.upsample+1+tau]];
                    [maxD maxClasses] = max(D,[],1);          
                    if self.P.storeSICvaluesForDebugging
                        SICvals.D(:,:,count) = D;
                        SICvals.subtractor(:,:,count) = subtractor;
                        SICvals.prior(count) = log(p(maxC));
                    end
                else
                    fprintf('Detection cancelled due to norm diff! (Y: %f - Y+sub: %f)\n', norm(Y), norm(Y+subtractor));
                    break
                end
            end
            if isempty(gdf)
                fprintf('No spikes in Epoch\n');
            end
        end
        
        %%% ------------------------------------------------------
        %%% -------------------PLOTTER----------------------------
        %%% ------------------------------------------------------
        function plotLastChunkSorting(self, varargin)
            import mysort.*
            P.start = [];
            P.stopp = [];
            P.titles = 1;
            P.srate = 1;
            P.figureTitle = [];
            P.figureName = 'plotLastChunkSorting';
            P.figureHandle = [];
            P.removeSpikesFromD = false;
            P.gui = false;
            P.gtGdf = [];
            P.X = [];
            P = util.parseInputs(P, 'plotLastChunkSorting', varargin);
            G = mysort.plot.globals();
            if isempty(P.start)
                P.start = self.chunk_start;
            end
            if isempty(P.stopp)
                P.stopp = self.chunk_end;
            end
            X = self.DH(P.start:P.stopp,:)';

            if P.removeSpikesFromD
                D = self.lastChunkSortingD;
            else
                D = self.D;
            end
            D = D(:, P.start-self.chunk_start_overlap+1 : P.stopp-self.chunk_start_overlap+1);
            
            % PLOT DATA WITH TEMPLATES SUPERIMPOSED
            if P.gui
                clf;
            else
                if isempty(P.figureHandle)
                    fig1 = mysort.plot.figure('color','w');
                else
                    fig1 = P.figureHandle;
                end
            end
            mysort.plot.figureName('BOTM.plotLastChunkSorting');
            
            if ~isempty(P.figureName)
                set(gcf, 'Name', P.figureName);
            end
            
            ax = mysort.plot.subplots([2,1], 'figureTitle', P.figureTitle,...
                            'offsetY',.08);
            axes(ax(1));
            % Plot data
            spacer = plot.mc(X, 'figure', 0, 'color', {'k'},...
                            'sampleOffset', P.start,...
                            'showZeroLine', true, 'srate', P.srate);
            hold on
            % Plot single spike templates
            for i=1:size(self.chunk_sorting,1)
                s1 = self.chunk_sorting(i,2);
                if s1>P.stopp ||  s1<P.start
                    continue
                end
                plot.mc(mysort.wf.v2m(self.T(self.chunk_sorting(i,1),:),size(self.DH,2)), 'figure', 0,...
                            'spacer',spacer,...
                            'linewidth', 2,...
                            'sampleOffset', self.chunk_sorting(i,2), ...
                            'color' , {plot.vectorColor(self.chunk_sorting(i,1))},...
                            'srate', P.srate);
            end
            % Plot single ground truth spikes if available
            if ~isempty(P.gtGdf)
                for i=1:size(P.gtGdf,1)
                    s1 = P.gtGdf(i,2);
                    if s1>P.stopp ||  s1<P.start
                        continue
                    end
                    plot.mc(mysort.wf.v2m(self.T(P.gtGdf(i,1),:),size(self.DH,2)), 'figure', 0,...
                                'spacer',spacer,...
                                'linewidth', 3, ...
                                'lineStyle', '--',...
                                'sampleOffset', P.gtGdf(i,2), ...
                                'color' , {plot.vectorColor(P.gtGdf(i,1))},...
                                'srate', P.srate);
                end                
            end
            axis tight
            
            % PLOT DISCRIMINATING FUNCTIONS
            axes(ax(2));
            hold on
            % Plot Log Prior for overlapping spikes
            for i=size(self.T,1)+1:size(self.D,1)
                plot.mc(D(i,:), 'figure', 0, 'spacer',0,...
                    'color', {'m'},'sampleOffset', P.start,...
                    'linewidth', 1, 'lineStyle', ':', 'srate', P.srate);
            end
            % Plot Log Prior of single spikes
            for i=1:size(self.T,1)
                plot.mc(D(i,:), 'figure', 0, 'spacer',0,...
                    'color', {plot.vectorColor(i)},...
                    'sampleOffset', P.start, 'linewidth', 2, 'srate', P.srate);
            end

            % Plot Noise Prior
            axis tight
            lims = get(gca, 'xlim');
            plot(lims, repmat(self.threshold,1,2), ':k');   
            
            linkaxes(ax, 'x'); 
            axis tight
            
            if P.titles
                title(ax(1), 'Data and superimposed templates', 'fontsize', G.title_fontsize);
                title(ax(2), 'Filter functions', 'fontsize', G.title_fontsize);
            end

            if P.gui
                fprintf(' Paused because of plot. Press button to continue.\n');
                pause
            end
        end
    end
    
    %%% ------------------------------------------------------
    %%% -----------------STATIC FUNCTIONS---------------------
    %%% ------------------------------------------------------    
    methods(Static)
        %%% ------------------------------------------------------        
        function logPr = calculateProbability(MAHA, priors)  
            %  probabilty calculation
            logPr = -.5*MAHA + repmat(log(priors),1,size(MAHA,2));
        end
        
        %%% ------------------------------------------------------        
        function post = calculatePosterior(logPr)          
            % Posterior calculation
            post = logPr./repmat(sum(logPr,1),size(logPr,1),1);            
        end        
    end % methods
end % classdef