
classdef EMClustering < mysort.sorters.GaussianNoiseSpikeSorterInterface
    properties (Constant=true)
        minNoiseEpochLength = 100;
    end
    properties
        spikesX
        pwSpikesX
        noiseSpikesX
        noiseMean
        spikesUnaligned
        fetX
        noiseFetX
        noiseTimes
    end
    
    methods
        %%% ------------------------------------------------------ 
        function self = EMClustering(varargin)
            self = self@mysort.sorters.GaussianNoiseSpikeSorterInterface(varargin{:});
        end
        
        %% ----------------CONSTRUCTOR INIT---------------------- 
        function init(self, varargin)
            init@mysort.sorters.GaussianNoiseSpikeSorterInterface(self, varargin{:});               
            self.P.detectArtefacts = 1;
            self.P.preWhiten    = 1;
            self.P.align        = 1;
            self.P.repeats      = 3;
            self.P.nClusterDims = 4;
            self.P.EXPSIG       = 6;
            self.P.nMinCluster  = 1;
            self.P.nMaxCluster  = 10;
            self.P.debug        = false;
            self.P = mysort.util.parseInputs(self.P, 'EMClustering', varargin, 1);
            
            self.spikesX = [];
            self.fetX = [];
            self.bReadyToSort = true;
        end
        
        %%% ------------------------------------------------------ 
        function sorting = sort_(self, X)
            [self.spikesX self.pwSpikesX times] = self.prepareSpikes(self.DH);
            
            [self.fetX IDs] = self.clusterSpikes(self.pwSpikesX);

            [self.spikesX, sorting, self.templates, self.fetX, self.noiseFetX,...
                self.noiseSpikesX, self.noiseTimes] =...
                    self.postClustering(self.spikesX, self.fetX, times, IDs);
        end
        
        %%% ------------------------------------------------------ 
        function [spikesX pwSpikesX times] = prepareSpikes(self, X)
            if self.bSorted
                spikesX = self.spikesX;
                pwSpikesX = self.pwSpikesX;
                times = self.sorting(:,1);
                C = self.NE.getNoiseCov();
                return
            end
            import mysort.*
            X = self.getX(X);
            spikeTimes = self.P.spikeDetector.detect(X); 
            spikeEpochs = [spikeTimes-self.P.cutleft spikeTimes-self.P.cutleft+self.P.Tf-1];
            spikeEpochs = epoch.removeOverlapping(spikeEpochs, self.artefactEpochs);

            % Remove invalid epochs
            nS_org = size(spikeEpochs,1);
            spikeEpochs = spikeEpochs(spikeEpochs(:,1)>0,:);
            spikeEpochs = spikeEpochs(spikeEpochs(:,2)<size(X,2),:);
            nS = size(spikeEpochs,1);

            if nS_org ~= nS
                fprintf('Warning! Invalid SpikeEpochs removed!! (extractSpikes, %d)\n', nS_org-nS);
            end    
            if nS == 0
                error('No Spikes were found in initialization!!!');
            elseif  nS < 10
                fprintf('Warning in stdSpikeSorting: Very few spikes found (%d)!\n', nS);
            else 
                fprintf('%d Spikes were found in initialzation.\n', nS);
            end    

            
            C = self.NE.getNoiseCov(X);
            iU = self.NE.iU;
            UP = 4;
            spikesX = mysort.epoch.extractWaveform(X, spikeEpochs);
            spikesX = mysort.util.resampleMC(spikesX, UP, 1);
            [tau spikesX] = mysort.util.alignWaveformsOnMaxOrMin(spikesX, self.nC);%, 'nIter', 10);  
%             tau = mysort.util.alignWaveformsOnMax(resPwSpikesX, self.nC);%, 'nIter', 10);
%             aliResPwSpikesX = mysort.util.shiftRows(resPwSpikesX, tau, self.nC, 1);           
            spikesX = mysort.util.resampleMC(spikesX, 1, UP);
            spikeEpochs = spikeEpochs - repmat(round(tau/UP),1,2);
            times = spikeEpochs(:,1);
            pwSpikesX = spikesX*iU;
            self.spikesX = spikesX;
            self.pwSpikesX = pwSpikesX;
        end
        
        %%% ------------------------------------------------------
        function [spikesX pwSpikesX times C] = getPreparedSpikes(self)
            assert(self.bSorted, 'Sorter did not prepare spikes yet!');
            spikesX = self.spikesX;
            pwSpikesX = self.pwSpikesX;
            times = self.sorting(:,2);
            C = self.NE.C;
        end
        
        %%% ------------------------------------------------------
        function [fetX IDs] = clusterSpikes(self, pwSpikesX)
            fetX = mysort.util.dimReductionPCA(pwSpikesX, self.P.nClusterDims);
            IDs = mysort.clustering.gmmClu(fetX,'kmin', self.P.nMinCluster,...
                                'kmax', self.P.nMaxCluster,....
                                'repeats', self.P.repeats,...
                                'EXPSIG', self.P.EXPSIG);
        end
        
        %%% ------------------------------------------------------
        function [spikesX, sorting, templates, fetX, noiseFetX, noiseSpikesX, noiseTimes] = ...
                postClustering(self, spikesX, fetX, times, IDs)
            % Remove Noise Spikes
            noiseID = 1;
            noiseSpikesX = spikesX(IDs==noiseID,:);
            spikesX = spikesX(IDs~=noiseID,:);
            
            noiseTimes = times(IDs==noiseID);
            times = times(IDs~=noiseID);
            
            noiseFetX = fetX(IDs==noiseID,:);
            fetX = fetX(IDs~=noiseID,:);
            IDs = IDs(IDs~=noiseID)-1;

            % Remove potential empty clusters:
            uClasses = unique(IDs);
            for i=1:length(IDs)
                IDs(i) = find(uClasses==IDs(i));
            end
            
            clus = unique(IDs);
            templates = zeros(length(clus), size(spikesX,2));
            for i=1:length(clus)
                templates(i,:) = mean(spikesX(IDs==clus(i),:),1);
            end
            sorting = [IDs times];         
        end       
        
        %%% ------------------------------------------------------ 
        function plotClustering(self)
            classes = self.getSpikeClasses();
            mysort.plot.clustering(self.fetX, classes);            
        end
        
        %%% ------------------------------------------------------ 
        function debugPlot(self)
            if self.P.debug
                mysort.plot.spikes(resPwSpikesX, 'nC', self.nC);
                title('Prewhitened spikes');
                mysort.plot.spikes(resPwSpikesX, 'nC', self.nC, 'classes', IDs);
                title('Prewhitened spikes, per cluster');
                mysort.plot.spikes(aliResPwSpikesX, 'nC', self.nC);
                title('Aligned prewhitened spikes');
                mysort.plot.spikes(aliResPwSpikesX, 'nC', self.nC, 'classes', IDs);
                title('Aligned prewhitened spikes, per cluster');                
            end
        end
    end
    
    methods (Static)

    end
end