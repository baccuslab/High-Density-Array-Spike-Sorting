classdef OfflineSpikeSorter < mysort.util.HandleObject
    properties (SetAccess=private)
        
    end
    properties
        P
        nC
        gdfs
        chunkCounter
        noise
        bufferSpikes
        bufferSpikeTimes
        bufferSpikesUpsampled
        bIsTrained
        
        featurePCs
        
        templates 
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = OfflineSpikeSorter(varargin)
            self.P.minDistFromSpikes = 60;
            self.P.maxLength = 200000;
            
            self.P.spikeCutting.cutLeft = 10;
            self.P.spikeCutting.Tf = 35;
            
            self.P.alignment.upsample = 4;
            self.P.alignment.upsample_cutEnds = 2;
            
            self.P.featureExtraction.Tf = 15;
            self.P.featureExtraction.nDims = 6;
            self.P.featureExtraction.cutLeft = 5;
            
            self.P.clustering.minSpikes = 100000;
            self.P.clustering.maxSpikes = 110000;
            self.P.clustering.minSpikesPerCluster = 10;
            self.P.clustering.bandwidthFactor = 1.5; 
            
            self.P = mysort.util.parseInputs(self.P, varargin, 'error');
            
            self.chunkCounter = 0;
            self.bIsTrained = 0;
            self.bufferSpikes = {};
            self.bufferSpikeTimes = {};
            self.bufferSpikesUpsampled = {};
        end
        
        % --------------------------------------------------------
        function processChunk(self, X, spikeTimes)
            
            % create handle class to avoid data copying
            X = mysort.ds.Matrix(X);            
            
            if isempty(self.noise)
                self.nC = size(X,2);
                self.learnNoiseCov(X, spikeTimes);
            end
            
            % cutspikes
            spikes_cut = X.getWaveform(spikeTimes,...
                         self.P.spikeCutting.cutLeft, self.P.spikeCutting.Tf);    
                     
            if self.bIsTrained
                self.processChunk_(spikeTimes, spikes_cut);
            else
                self.bufferChunk(X, spikeTimes, spikes_cut);
            end
        end
        
        % --------------------------------------------------------
        function processChunk_(self, spikeTimes, spikes_cut)
            self.chunkCounter = self.chunkCounter+1;
            disp('Classifying');
            S = spikes_cut;
            T = self.templates.wfs;
            D = bsxfun(@plus, full(dot(S',S',1))', full(dot(T',T',1))) - full(2*(S*T'));
            [m ids] = min(D, [], 2);      
            
%             [D maxTaus TupShiftDown FupShiftDown tauRange tauIdx subTauRange subTauIdx] = mysort.wf.vTemplateMatching(...
%                     S.spikeCut.wfs, T, nC, idxstart, 'maxShift', 5, 'upsample', 5, 'noiseCovariance', C);
%                 F = FupShiftDown(:,:,1);
%                 EF = diag(T*F');                     % compute energies
%                 Prior  = P.templateMatchingCut.prior;                
%                 DISCR = D - .5 * repmat(EF', nSpCut, 1)  + log(Prior);  % compute botm dicriminant            
            
            chunk_gdf = [ids(:) spikeTimes(:)];
            self.gdfs{self.chunkCounter} = chunk_gdf;
        end
        
        % --------------------------------------------------------
        function train(self)
            disp('Training');
        	vS = cell2mat(self.bufferSpikesUpsampled(:));
            tS = mysort.wf.v2t(vS, self.nC);
        
            [tS tau_m mMPmask] = mysort.wf.tAlignOnMean(tS,...
                'restrictToNMaximalValues', 10*self.P.alignment.upsample*self.nC,...
                'meanMaskMinLen', 8*self.P.alignment.upsample,...
                'maxShift', 3*self.P.alignment.upsample,...
                'maxShiftPerIter', self.P.alignment.upsample*2,...
                'useMedian', true, ...
                'maxIter', 5);

            % Cut ends away again and downsample
            selIdx = (self.P.spikeCutting.cutLeft-self.P.featureExtraction.cutLeft+1-self.P.alignment.upsample_cutEnds);
            selIdx = (selIdx:selIdx+self.P.featureExtraction.Tf-1)*self.P.alignment.upsample;
            vtS = mysort.wf.t2v(tS);
            tS = tS(selIdx,:,:);
            
            % Prewhiten
            
            Spw = mysort.wf.t2v(tS)/self.noise.U;
            
            % PCA
            [F self.featurePCs] = mysort.util.dimReductionPCA(Spw,...
                    self.P.featureExtraction.nDims, [], 3*1000000);
                
            % clustering
            bw = sqrt(self.P.clustering.bandwidthFactor * self.P.featureExtraction.nDims);
            [clustCent, point2cluster, clustMembsCell] = MeanShiftCluster(F', bw);
            [ids_clu T] = MeanShiftClusterBundleResult(F, clustMembsCell, self.P.clustering.minSpikesPerCluster);
            
            self.templates.wfs_ali = mysort.util.calculateClassMeans(vtS, ids_clu, 'usemedian');
            self.templates.wfs_pw = mysort.util.calculateClassMeans(Spw, ids_clu, 'usemedian');
            self.templates.wfs_up = mysort.util.calculateClassMeans(vS, ids_clu, 'usemedian');
            clear vtS Spw vS tS 
            vS = cell2mat(self.bufferSpikes(:));
            self.templates.wfs = mysort.util.calculateClassMeans(vS(1:length(ids_clu),:), ids_clu, 'usemedian');
            self.templates.ids = unique(ids_clu);
            
            self.bufferSpikesUpsampled = {};
            
            for i=1:length(self.bufferSpikes)
                self.processChunk_(self.bufferSpikeTimes{i}, self.bufferSpikes{i});
            end
            self.bufferSpikes = {};
            self.bufferSpikeTimes = {};
            self.bIsTrained = true;
        end
        
        % --------------------------------------------------------
        function bufferChunk(self, X, spikeTimes, spikes_cut)
                     
            % buffer all cut spikes
            self.bufferSpikes{end+1} = spikes_cut;
            self.bufferSpikeTimes{end+1} = spikeTimes;             

            % upsample only so many spikes as needed
            nSpikesInBuffer = sum(cellfun(@(x) size(x,1), self.bufferSpikes));
            
            % only process as many spikes as needed for clustering
            if nSpikesInBuffer > self.P.clustering.maxSpikes
                nExcessSpikes = nSpikesInBuffer-self.P.clustering.maxSpikes;
                keepSpikes = min(size(spikes_cut,1),max(10,size(spikes_cut,1)-nExcessSpikes+1));
                spikes_cut = spikes_cut(1:keepSpikes,:);
            end            
            spikes_cut = mysort.wf.tResample(mysort.wf.v2t(spikes_cut, size(X,2)), self.P.alignment.upsample, 1);
            
            % remove artefacts not needed for alignment
            spikes_cut([1:self.P.alignment.upsample_cutEnds*self.P.alignment.upsample ...
                        end-self.P.alignment.upsample_cutEnds*self.P.alignment.upsample+1:end],:,:) = [];
            spikes_cut = mysort.wf.t2v(spikes_cut);
            
            self.bufferSpikesUpsampled{end+1} = spikes_cut;
            
            if nSpikesInBuffer >= self.P.clustering.minSpikes
                self.train();
            end
        end
        
        % --------------------------------------------------------
        function learnNoiseCov(self, X, spikeTimes)
            % estimate noise cov
            s1 = spikeTimes - self.P.minDistFromSpikes;
            s2 = spikeTimes + self.P.minDistFromSpikes;
            L  = min(self.P.maxLength, size(X,1));
            epochs = mysort.epoch.flip(mysort.epoch.merge([s1 s2]), L);

            Cest = mysort.noise.Covest2(X, 'maxLag', self.P.featureExtraction.Tf,...
                'maxSamples', L, 'noiseEpochs', epochs, 'forceMethod', 'xcorr');
            
            self.noise.C_unloaded = mysort.noise.ccol2Cte(Cest.CCol, self.P.featureExtraction.Tf);
            self.noise.meanNoiseStd = sqrt(mean(diag(self.noise.C_unloaded)));
            self.noise.CestS = Cest.toStruct();  
            self.noise.ccol = self.noise.CestS.CCol;
            
            % Prepare Noise Covariance matrix
            nC = size(X,2);
            self.noise.ccol(1:nC, 1:nC) = self.noise.ccol(1:nC, 1:nC) + ...
                                           diag(diag(self.noise.ccol(1:nC, 1:nC)));
            self.noise.ccol = self.noise.ccol/2;
            self.noise.C = mysort.noise.ccol2Cte(self.noise.ccol, self.P.featureExtraction.Tf);
            self.noise.U = chol(self.noise.C);        
        end
    end
end
          

    
    