classdef BOTMOvp < mysort.sorters.FilterBasedSpikeSorterInterface
    
    properties (SetAccess=private)
        init_counter
    end
    properties
        E
        D
        DOvpC
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
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = BOTMOvp(Covest, T, varargin)
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
            self.P = mysort.util.parseInputs(self.P, varargin);
            
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
%             caller_name = mysort.util.getCallerName()
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
             [sorting filterOut] = self.sortData(t1, t2);
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
            self.DOvpC = mysort.sorters.DiscriminantFunctionContainer(self.D', self.CONF, 6);
            fout = [self.D; self.DOvpC.DOvp'];
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
            [self.T self.templateAlignmentShifts] = mysort.wf.tAlignOnCorrelation(tT, 'trunc', 1, 'absMax', 1);    
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
            self.CONF_up = mysort.util.resampleTensor(self.CONF, self.P.upsample, 1);

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
            
            % For an Overlap D-Fun to win, it must be higher then BOTH
            % individual D-Funs at the corresponding timeshifts!
            [mPr c] = max([self.D; self.DOvpC.DOvp'],[],1);
            [amps, idc] = findpeaks(mPr,'MINPEAKHEIGHT',self.threshold, 'MINPEAKDISTANCE', 8);
            
            % Check if an overlap function was detected. We only accept
            % overlap detections if at the corresponding points the single
            % spike discriminant functions are lower than the overlap
            % function. Otherwise a spike of B might be detected as an
            % overlap between A and B with shift t because it is more
            % likely to have A+B rather than A since there is no A but B.
            % If this function peaks with distance t from the peak of
            % function B it might be detected as a spike
            chunk_sorting = [c(idc(:))' idc(:)];
            ovpDetections = find(chunk_sorting(:,1)>size(self.D,1));
            for i=length(ovpDetections):-1:1
                ovpid   = c(idc(ovpDetections(i)));
                ovptime = idc(ovpDetections(i));
                ovpamp  = amps(ovpDetections(i));
                ovpidxidx = self.DOvpC.DOvpIndex(:,1) == ovpid;
                ovptype = self.DOvpC.DOvpIndex(ovpidxidx,:);
                tIdx = ovptime+ovptype(4);
                if size(self.D,2) < tIdx || ...
                   tIdx < 1 || ...
                   self.D(ovptype(2), ovptime) >= ovpamp || ...
                   self.D(ovptype(3), tIdx) >= ovpamp
                    chunk_sorting(ovpDetections(i),:) = [];
                end
            end
            if 0
                figure;
                ah = subplot(2,1,1);
                plot(self.D');
                ah(2) = subplot(2,1,2);
                plot(self.DOvpC.DOvp);
                linkaxes(ah);                
            end
            if ~isempty(chunk_sorting)
                chunk_sorting = sortrows(chunk_sorting,2);
            end
            chunk_sorting(:,2) = chunk_sorting(:,2)-ceil(self.Tf/2)-1;
            [chunk_sorting wasResolved] = mysort.spiketrain.resolveOvpIndexInGdf(chunk_sorting, self.DOvpC.DOvpIndex);
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