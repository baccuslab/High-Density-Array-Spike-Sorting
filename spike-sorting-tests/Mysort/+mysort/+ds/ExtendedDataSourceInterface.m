classdef ExtendedDataSourceInterface < mysort.ds.DataSourceInterface
    properties
        memoryBufferNoiseSmad
        spikeTrains
    end
    
    methods(Abstract)
        getData_(self, idx1, idx2)
    end
    
    methods
        %------------------------------------------------------------------
        function self = ExtendedDataSourceInterface(varargin)
            self = self@mysort.ds.DataSourceInterface(varargin{:});
            self.memoryBufferNoiseSmad = [];
        end
        
        %------------------------------------------------------------------
        function Cest = getCovest(self, maxlag, maxsamples, maxdist, forceMethod)
            if nargin < 2 
                maxlag = 79;
            end
            if nargin < 3
                maxsamples = 150000;
            end
            if nargin < 4
                maxdist = 40;
            end
            if nargin < 5
                forceMethod = []; %'matmul', 'xcorr'
            end
            fprintf('Calculating Covest, that may take a while...\n');
            [times pks] = self.detectSpikes();
            times = double(cell2mat(times))';
            spikeEpochs = mysort.epoch.merge([times(:)-50 times(:)+50]);
            noiseEpochs = mysort.epoch.flip(spikeEpochs, size(self,1));
            t1 = tic;
            Cest = mysort.noise.Covest2(self, 'maxLag', maxlag, ...
                'maxSamples', maxsamples, 'noiseEpochs', noiseEpochs,...
                'maxDist', maxdist, 'forceMethod', forceMethod);
            t2 = toc(t1);
            disp('Done.'); disp(t2);
        end
        %------------------------------------------------------------------
        function R = xcorr(self, varargin)
            P.channelIdx = 1:self.size(2);
            P.maxLen = 200000;
            P.maxLag = 100;            
            P.normalization = 'none';
            P = mysort.util.parseInputs(P, varargin);
            
            R = xcorr(self(1:P.maxLen, P.channelIdx), P.maxLag, P.normalization);            
        end
        
        %------------------------------------------------------------------
        function [smad] = noiseStd(self, varargin)
            % Calculate channel wise noise standard deviation with the median
            % absolute deviation (MAD), invert data to ignore negative peaks
            % for that calculation
            P.channelIdx = 1:self.size(2);
            P.maxLen = 300000;
            P.thr = 4;
            P.Tf = 80;            
            P = mysort.util.parseInputs(P, varargin);
           
            Len = self.size(1);
            fullChanIdx = P.channelIdx;
            if isempty(self.fullMultiElectrode)
                nC = self.MultiElectrode.getNElectrodes();
            else
                nC = self.fullMultiElectrode.getNElectrodes();
                if ~isempty(self.activeChannels)
                    fullChanIdx = self.activeChannels(P.channelIdx);
                end
            end
            
            if isempty(self.memoryBufferNoiseSmad)
                self.memoryBufferNoiseSmad = nan(1, nC);
            end
            notCalcIdx = isnan(self.memoryBufferNoiseSmad(fullChanIdx));
            if any(notCalcIdx)
                cidx = P.channelIdx(notCalcIdx);
                fullcidx = fullChanIdx(cidx);
                disp('Computing noise std...'); tic
                smadL = min(Len, P.maxLen);
                smad = mysort.noise.estimateSigma(...
                        self.getData(1:smadL, cidx), P.Tf, P.thr);
                self.memoryBufferNoiseSmad(fullcidx) = smad;
                disp('Done.'); toc        
            end
            smad = self.memoryBufferNoiseSmad(fullChanIdx);
        end
        %------------------------------------------------------------------
        function [times pks] = detectSpikes(self, varargin)
            P.channelIdx = 1:self.size(2);
            P.chunkSize = 100000;
            P.thr = 3.5;
            P.energyfun = @(x) -x;
            P.minPeakDistance = ceil(self.getSamplesPerSecond/1000); % 1ms
            P.Len = [];
            P = mysort.util.parseInputs(P, varargin);
            
            if isempty(P.Len)
                P.Len = self.size(1);
            end
            
            % get noise std
            smad = self.noiseStd('channelIdx', P.channelIdx);
            
            % Detect spikes in the beginning of the file
            disp('Detecting spikes...'); tic
            pks = cell(length(P.channelIdx),1);
            times = pks;
            for cidx = 1:length(P.channelIdx)
                c = P.channelIdx(cidx);
                pks{c,1} = [];
                times{c,1} = [];
            end
            chunker = mysort.util.Chunker(P.Len, 'chunkSize', P.chunkSize, ...
                'progressDisplay', 'console', 'minChunkSize', 1000, 'chunkOverlap', 2*P.minPeakDistance);
            while chunker.hasNextChunk()
                [chunkOvp chunk] = chunker.getNextChunk();
                X = double(self.getData(chunkOvp(1):chunkOvp(2), P.channelIdx));
                for cidx = 1:length(P.channelIdx)
                    c = P.channelIdx(cidx);
                    [pks_, times_] = findpeaks(P.energyfun(X(:,c)), 'MINPEAKHEIGHT', smad(cidx)*P.thr,...
                                                            'MINPEAKDISTANCE', P.minPeakDistance);
                    pks_ = X(times_,c); % get the right amplitudes! (sign!)
                    pks_ = pks_(:);
                    times_ = times_(:);
                    % remove spikes that are outside this chunk
                    rmvIdx = (times_+chunkOvp(1) < chunk(1)) | (times_+chunkOvp(1) > chunk(2));
                    pks_(rmvIdx) = [];
                    times_(rmvIdx) = [];
                    
                    pks{c,1} = [pks{c}; pks_];
                    times{c,1} = [times{c}; times_+chunkOvp(1)-1];
                end
            end
            disp('Done.'); toc    
        end
        %------------------------------------------------------------------
        function allspikes = getMergedSingleElectrodeDetectedSpikes(self, mergeSpikesMaxDist, varargin)
            [times pks] = self.detectSpikes(varargin{:});
            allspikes = sortrows([cell2mat(times)   cell2mat(pks)], 1);
            allspikes  = mysort.spiketrain.mergeSingleElectrodeDetectedSpikes(allspikes, mergeSpikesMaxDist);
        end        
    end
end

%     P.cut_Tf = 30;
%     P.cut_left = 5;
%     P.cut_upsample = 8;
    
%     % Cut spikes and upsample
%     disp('Cutting and upsampling spikes...'); tic
%     unlimited = H5ML.get_constant_value('H5S_UNLIMITED');
%     dims = [1 P.cut_Tf*P.cut_upsample];
%     maxDims = [unlimited P.cut_Tf*P.cut_upsample];
%     h5type = 'H5T_NATIVE_DOUBLE';
%     chunk_dims = dims;
%     deflation = 0;
%     spike_fname = [P.outFileName '_spks.h5'];
%     spike_h5path = '/cutspikes';
%     SP = mysort.h5.matrixCreate(spike_fname, spike_h5path, dims, maxDims,...
%         h5type, chunk_dims, deflation);
%     chan_time_amps_h5path = '/chan_time_amps';
%     CHAN_TIME_AMPS = mysort.h5.matrixCreate(spike_fname, chan_time_amps_h5path,...
%         [1 5], [unlimited 5], h5type, [1000 5], deflation);
%     sidx = 1;
%     for c=1:nC
%         for s=1:length(mi_pks{c,1})
%             s1 = max(mi_pks{c,1}(s)-P.cut_left, 1);
%             s2 = min(mi_pks{c,1}(s)-P.cut_left+P.cut_Tf-1, size(MEA,1));
%             sp = double(MEA(s1:s2,c));
%             rsp = resample(sp, P.cut_upsample, 1);
%             [mi mi_idx] = min(rsp);
%             [ma ma_idx] = max(rsp);
%             SP(sidx, 1:length(rsp)) = rsp';
%             real_mi_time = s1+(mi_idx-1)/P.cut_upsample;
%             real_ma_time = s1+(ma_idx-1)/P.cut_upsample;
%             mi_pks{c,1}(s) = mi;
%             ma_pks{c,1}(s) = ma;
%             times{c,1}(s) = real_mi_time;
%             CHAN_TIME_AMPS(sidx, :) = [c real_mi_time mi real_ma_time ma]; 
%             sidx=sidx+1;
%         end
%     end
%     disp('Done.'); toc

%     % Build sorted matrix containing time points, channel index and
%     % amplitude of found spikes
%     disp('Building spike index...'); tic
%     times(cellfun(@isempty, times)) = [];
%     mi_pks(cellfun(@isempty, mi_pks)) = [];
%     allTimes = cell2mat(times);
%     allAmps   = cell2mat(mi_pks);
%     allChans  = zeros(size(allTimes));
%     nCum = 0;
%     for i=1:length(times)
%         allChans(nCum+1:nCum+length(times{i}),1) = i;
%         nCum = nCum + length(times{i});
%     end
%     timesChansAmpsSorted = sortrows([allTimes allChans allAmps]);
%     fprintf('maxtime: %d\n', max(allTimes));
%     disp('Done.'); toc