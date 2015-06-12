
classdef BOTMWA < mysort.sorters.BOTM
    
    properties (SetAccess=private)
        
    end
    properties
        D_
        optimal_alphas
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = BOTMWA(varargin)
            self = self@mysort.sorters.BOTM(varargin{:});
        end
        
       
        
        %%% ------------------------------------------------------
        function fout = getFilterOutputs_(self, X, filterindex)
            self.debugout('Applying filter...', self.LEVEL_PROCESSING);
            self.Y = mysort.util.applyMcFilters(X, self.F);

            self.debugout('Calculating Discriminant Function...', self.LEVEL_PROCESSING);
%             [self.D self.Dshifts] = mysort.util.calculateDiscriminantFunctions(self.Y, self.CONF, self.priors);  
            [self.D_ self.optimal_alphas] = mysort.util.calculateDiscriminantFunctionsWithAlpha(self.Y, self.CONF, self.priors);  
            Palpha = evpdf(self.optimal_alphas, 1, .2);
            Palpha(self.optimal_alphas<.1) = 0;
            logPalpha = log(Palpha);
            logPalpha(logPalpha < -2) = nan;                
            self.D = self.D_+logPalpha;
            if 1
                [D Dshifts] = mysort.util.calculateDiscriminantFunctions(self.Y, self.CONF, self.priors); 
                ax = [];
                mysort.plot.figure('w', 800, 'h', 800);
                ax(1) = subplot(4,1,1)
                hold on
                for i=1:size(self.D,1)
                    plot(D(i,:), 'color', mysort.plot.vectorColor(i), 'linewidth', 2);
                end
                title('BOTM discrm func')
                

                ax(2) = subplot(5,1,2);
                hold on
                for i=1:size(self.D,1)
                    plot(self.optimal_alphas(i,:), 'color', mysort.plot.vectorColor(i), 'linewidth', 2);
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
                    plot(self.D_(i,:), 'color', mysort.plot.vectorColor(i), 'linewidth', 2);
                end
                title('Discr wo alpha prior')
                
                ax(5) = subplot(5,1,5)
                hold on
                for i=1:size(self.D,1)
                    plot(self.D(i,:), 'color', mysort.plot.vectorColor(i), 'linewidth', 2);
                end
                linkaxes(ax, 'x');
                [maxis maxChan] = max(self.D,[], 1);
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
        function chunk_sorting = Postprocessing(self)
            % Does the postprocessing of the filter calculation for this
            % CHUNK. D are the Discriminant functions of the current chunk,
            % CHUNK is the original data with overlapping regions to the
            % left and right CHUNK. 
            
            [mPr c] = max(self.D,[],1);
            spikeEpochs = mysort.epoch.fromBinaryVectorMinLen(...
                                    mPr>self.threshold, 150);%round(.5*self.Tf)
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
                    Yup(1:self.nF,:) = mysort.util.resampleMC(...
                        self.Y(1:self.nF,epochs(e,1)-o1:epochs(e,2)+o2), self.P.upsample, 1);
                    Yup = Yup(:, o1*self.P.upsample+1:end-o2*self.P.upsample);
%                     Aup(1:self.nF,:) = mysort.util.resampleMC(...
%                         self.optimal_alphas(1:self.nF,epochs(e,1)-o1:epochs(e,2)+o2), self.P.upsample, 1);
%                     Aup = Aup(:, o1*self.P.upsample+1:end-o2*self.P.upsample);
                else
                    Yup = self.Y(:,epochs(e,1):epochs(e,2));
%                     Aup = self.optimal_alphas(:,epochs(e,1):epochs(e,2));
                end
                epoch_gdf = self.processEpoch(Yup, self.CONF_up, self.Tf_up, self.priors, epochs(e,:));
                if ~isempty(epoch_gdf)
                    epochOffset = epochs(e,1)-ceil(self.Tf/2)-1;
                    epoch_gdf(:,2) = round(epoch_gdf(:,2)./self.P.upsample) + epochOffset;
                    chunk_gdf = [chunk_gdf; epoch_gdf];
                end
            end
        end

        %%% ------------------------------------------------------   
        function gdf = processEpoch(self, Y, M, Tf, p, epoch)
            % Detect and remove Spikes as long as discrimant functions
            % are larger than the noise prior
            % - uses SIC
%             D = mysort.util.calculateDiscriminantFunctions(Y, M, p);
            gdf = [];
            count = 0;
            debug = 1;
            if debug
%                 close all
                fig = mysort.plot.figure('h', 1000);
                ah = mysort.plot.subplot([4,1]);
                        
                %mysort.plot.XIvsF(randn(3,Tf-2), randn(3,Tf-2), 'XIvsF', self.CONF_up)
            end
            subtractor = [];
%             D_old = [];
            while 1
                [D A] = mysort.util.calculateDiscriminantFunctionsWithAlpha(Y, M, p);  
                Palpha = evpdf(A, 1, .2);
                Palpha(A<.1) = 0;
                logPalpha = log(Palpha);
                logPalpha(logPalpha < -2) = nan;                
                D = D+logPalpha;                
                [maxD maxClasses] = max(D,[],1); 
                if ~any(maxD>self.threshold)
                    break
                end
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
                    sub = A(maxC, t) * squeeze(M(:,f,maxC))';
%                     if f == maxC
                        % Block double self detections?
%                         subtractor(f,:) = -inf;
%                     else                    
                        subtractor(f,:) = mysort.util.shiftSubtract(zeros(1,size(D,2)),...
                                        sub, offset, false);
%                     end
        
                end
                
                % Check if the subtraction of the expected response is
                % indeed decreasing the norm of the filter outputs.
                % Does not work for MEA data yes. Set "norm_constraint" to
                % "false"
                if ~self.P.norm_constraint || ...
                        norm(Y) > norm(Y+subtractor)
                    if debug
                        axes(ah(1));
                        title('X');
                        Margin = 200*self.P.upsample;
                        Margin_up = 200*self.P.upsample;
                        Tf_up = Tf*self.P.upsample;
                        
%                         xrange = Margin_up:self.P.upsample:epoch(2)-epoch(1)+Margin_up;
                        x = self.DH(max(1,epoch(1)-Margin_up):min(end, epoch(2)+Margin_up),:)';
                        plot(x,...
                            'k', 'linewidth', 3);
                        title('data (x)')
                        
                        axes(ah(2));
                        title(['D' num2str(count)]);
                        plot([1 size(D,2)], [0 0], ':k');
                        hold on                                        
                        if ~isempty(subtractor_old)
                            plot(D_old','--', 'linewidth',3);
                        end   
                        plot(D','linewidth',3);
                        title('D and old D (striped)')
%                         spacer = mysort.plot.mc(D,'plotZeroLine',1,...
%                             'linewidth',3,'figure',0,'color',{'k'}, 'figure', 0);
        
                    end
                    D_old = D;
                    Y = Y + subtractor; % + log(p(maxC));
                    if debug
                        %mysort.plot.mc(D, 'spacer', spacer,'linewidth',2 , 'color', {'r'},'figure',0);
                        plot(D',':', 'linewidth',2);
                        axes(ah(3));                        
                        title('Subtractor');
                        %mysort.plot.mc(subtractor + log(self.priors(maxC)), 'spacer', spacer,'linewidth',2,'figure',0,'color', {'g'});
                        plot((subtractor + log(self.priors(maxC)))', 'linewidth',2);
                        title('subtractor + priors')
                        
                        [newD newA] = mysort.util.calculateDiscriminantFunctionsWithAlpha(Y, M, p); 
                        Palpha = evpdf(newA, 1, .2);
%                         Palpha(newA<.1) = 0;
                        logPalpha = log(Palpha);
%                         logPalpha(logPalpha < -2) = nan;                
                        newD = newD+logPalpha;                            
                        axes(ah(4));                        
                        title('New D');
                        %mysort.plot.mc(subtractor + log(self.priors(maxC)), 'spacer', spacer,'linewidth',2,'figure',0,'color', {'g'});
                        plot(newD', 'linewidth',2);
                        linkaxes(ah, 'x');
                        axis tight
                    end
                    % correct for any shift that was done to the waveform
                    % in the initialization
                    tau = self.P.upsample*self.templateAlignmentShifts(maxC);
                    gdf = [gdf; [maxC t+self.P.upsample+1+tau A(maxC, t)]];                  
                else
                    fprintf('Detection cancelled due to norm diff! (Y: %f - Y+sub: %f)\n', norm(Y), norm(Y+subtractor));
                    break
                end
            end
            if isempty(gdf)
                fprintf('No spikes in Epoch\n');
            end
        end
    end
end % classdef