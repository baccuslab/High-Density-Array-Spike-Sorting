
classdef DetectorInterface < mysort.util.CacheableClass
    properties
        energy
        threshold
        times
        epochs
        DH
    end
    
    methods (Abstract)
        energy = computeEnergy(self, X)
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------     
        function self = DetectorInterface(varargin)
            self = self@mysort.util.CacheableClass(varargin{:});
        end
        
        %%% ----------------CONSTRUCTOR INIT---------------------- 
        function init(self, varargin)
            self.P.maximaWidth = 15;
            self.P.smoothingRange = 15;
            self.P.refractoryPeriod = 30;
            self.P.verbose = false;
            self.P.useMedian = false;
            self.P.threshFactor = 4;
            self.P.useThreshFactor = true;  
            self.P.degree = 2;
            self.P.threshold = [];
            suppress_warning = 1;
            self.P = mysort.util.parseInputs(self.P, 'DetectorInterface', varargin, suppress_warning);
            
            self.energy = [];
            self.threshold = [];
            self.times = [];
            %self.epochs = [];
        end
        
      
        %%% ------------------------------------------------------   
        function times = detect(self, X, srate)           
            if isa(X, 'mysort.datafile.DataFileInterface')
                self.DH = X;
            else
                if nargin< 3
                    srate = 32000;
                end                 
                self.DH = mysort.datafile.DataHandle(X, srate);
            end
            % get the number of channels and make sure X has the right format

            self.energy = self.computeEnergy(self.DH.getData());
            % do some smoothing and then thresholding to get out the spike times 
            
            [self.times self.threshold] = self.processEnergy(self.energy);  
            times = self.times;
        end    
        %%% ------------------------------------------------------   
        function epochs = detectEpochs(self, X, srate)
            if isa(X, 'mysort.datafile.DataFileInterface')
                self.DH = X;
            else
                if nargin< 3
                    srate = 32000;
                end 
                self.DH = mysort.datafile.DataHandle(X, srate);
            end
            % get the number of channels and make sure X has the right format

            self.energy = self.computeEnergy(self.DH.getData());
            % do some smoothing and then thresholding to get out the spike times 
            
            [self.epochs self.threshold] = self.processEnergy(self.energy, true);  
            epochs = self.epochs;
        end 
        %%% ------------------------------------------------------   
        function [times threshold] = processEnergy(self, energy, detectEpochs)
            % This function takes the energy of a detection method and
            % extracts occurence times from it. The output signal from the detector
            % may be one row vector or a matrix with outputs from each channel, i.e. a
            % [nrChannel x time] matrix.
            if nargin < 3
                detectEpochs = false;
            end
            nC = size(energy, 1);
            % threshold for the filter output, may either be a multiple of the output's
            % std or a fixed constant

            if self.P.useThreshFactor
                % compute a threshold based on the standard deviation of the spike
                % detection output signal
%                 if nC > 1
                    if isscalar(self.P.threshFactor)
                        self.P.threshFactor = ones(nC,1)*self.P.threshFactor;
                    elseif size(self.P.threshFactor,2) > size(self.P.threshFactor,1)
                        self.P.threshFactor = self.P.threshFactor'; % now a column vector
                    end
                    if self.P.useMedian 
                        threshold = median(energy,2);   % median of the channels in signal
                    else
                        threshold = std(energy,[],2);   % standard deviation of the channels in signal
                    end
                    threshold = threshold.*self.P.threshFactor; % + mean(signal, 2);
%                 else
%                     if self.P.useMedian
%                         threshold = median(energy)*self.P.threshFactor + mean(energy);
%                     else
%                         threshold = std(energy)*self.P.threshFactor + mean(energy);
%                     end
%                 end
            else
%                 if nC > 1
                    if isscalar(self.P.threshold)
                        threshold = self.P.threshold*ones(nC,1);
                    elseif size(self.P.threshold,2) > size(self.P.threshold,1)
                        threshold = self.P.threshold'; % now a column vector
                    else
                        threshold = self.P.threshold;
                    end
%                 else
%                     threshold = self.P.threshold;
%                 end
            end

            if ~detectEpochs
                % extract the spike time indices by thresholding the smoothed response
                times = mysort.detectors.DetectorInterface.applyThreshold(...
                    energy, threshold, self.P.maximaWidth);

                % remove spikes that are too close to their predecessor
                times = self.applyRefractoryPeriod(...
                    times, self.P.refractoryPeriod);
            else
                times = self.applyThresholdForEpochs(energy, threshold);
            end
        end
       
        
        %%% ------------------------------------------------------  
        function plotDetection(self, varargin)
            P.start = [];
            P.stopp = [];
            P = mysort.util.parseInputs(P, '', varargin);

            [X P.start P.stopp] = self.DH.getData('start', P.start, 'stopp', P.stopp);
            spacer = mysort.plot.mc(X,'sampleOffset', P.start, 'color', {'k'});
            hold on
            t = self.times(self.times>=P.start & self.times<=P.stopp);
            for i=1:size(X,1)
                s = (i-1)*spacer;
                plot(t, s+zeros(1,length(t)), '.g', 'markersize', 19);
                if ~isempty(self.threshold)
                    if strcmp(self.P.method, 'none') ||...
                       strcmp(self.P.method, 'abs') 
                        thr = self.threshold(i);
                        plot([P.start P.stopp], thr+repmat(s ,1,2), ':r');
                    end
                    if strcmp(self.P.method, 'minus') ||...
                       strcmp(self.P.method, 'abs') 
                        thr = -self.threshold(i);
                        plot([P.start P.stopp], thr+repmat(s ,1,2), ':r');
                    end 
                end
            end
            set(gca, 'xlim', [P.start P.stopp]);
            axis tight
        end
    end
    
    methods (Static)
        %%% ------------------------------------------------------  
        function epochs = applyThresholdForEpochs(energy, thresholds)
            if size(energy,1)>1 && size(thresholds,2)==1 && size(energy,2)>1
                thresholds = repmat(thresholds, 1, size(energy,2));
            end
            epochs = mysort.epoch.fromBinaryVector(sum(energy>thresholds)>0);
        end

        %%% ------------------------------------------------------
        function timeIndices = applyThreshold(energy, thresholds, windowlength)

            maximaIndices = [];
            for chan = 1:size(energy,1)
                [amps idx] = findpeaks(energy(chan,:),...
                    'minPeakHeight', thresholds(chan), 'minPeakDist', windowlength);
                maximaIndices = [maximaIndices idx];
            end
            maximaIndices = sort(maximaIndices);

            % remove spikes which are too near in more than one channel
            if length(maximaIndices)>1
                diffs = diff(maximaIndices);
                maximaIndices = [maximaIndices(1) maximaIndices([0 diffs] > windowlength)]';
            end
            timeIndices = maximaIndices;        
        end
        
        %%% ------------------------------------------------------  
        function spikeTimes = applyRefractoryPeriod(spikeTimes, refractoryPeriod)
            % Removes spikes that are too close to their predecessor. Computes the
            % difference between two consecutive spikes and deletes the second spike if
            % the time difference is smaller than the given refractoryPeriod. The
            % spikeTimes and the refractoryPeriod must be in the same units.

            i=1;
            while i<length(spikeTimes)
                if spikeTimes(i+1)-spikeTimes(i) < refractoryPeriod
                    spikeTimes(i+1) = [];
                else
                    i = i+1;
                end
            end
        end
    end
end