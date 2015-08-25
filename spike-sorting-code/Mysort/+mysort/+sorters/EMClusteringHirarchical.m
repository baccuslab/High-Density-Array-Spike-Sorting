
classdef EMClusteringHirarchical < mysort.sorters.EMClustering
    properties (Constant=true)
        
    end
    properties
        
    end
    
    methods
        %%% ------------------------------------------------------ 
        function self = EMClusteringHirarchical(varargin)
            self = self@mysort.sorters.EMClustering(varargin{:});
        end
        
        %%% ----------------CONSTRUCTOR INIT---------------------- 
        function init(self, varargin)
            init@mysort.sorters.EMClustering(self, varargin{:});              
            self.P.nSubDimensions = 2;
            self.P.initialSpikeSorter = [];
%             self.P.EXPSIG = 3;
            self.P = mysort.util.parseInputs(self.P, 'EMClusteringHirarchical', varargin, 1);
        end
        
        %%% ------------------------------------------------------ 
        function sorting = sort_(self, X)
            if isa(self.P.initialSpikeSorter, 'mysort.sorters.EMClustering')
                [spikesX pwSpikesX times self.C] = self.P.initialSpikeSorter.getPreparedSpikes();
                EM = self.P.initialSpikeSorter;
            else
                [spikesX pwSpikesX times self.C] = self.prepareSpikes(X);
                EM = sorters.EMClustering('otherP', self.P);
            end
            EM.sort();            
            initClasses  = EM.getSpikeClasses();         
            
            IDs = self.reClusterSpikes(pwSpikesX, initClasses);
            
            [self.spikesX, sorting, self.templates, self.fetX, self.noiseFetX,...
                self.noiseSpikesX, self.noiseTimes] =...
                    self.postClustering(spikesX, EM.fetX, times, IDs);
        end
        %%% ------------------------------------------------------
        function [classes] = reClusterSpikes(self, pwSpikesX, classes)
            import mysort.*

            cls = unique(classes);
            nClasses = length(cls);
            % Compute candidate clustering projections as those dimensions
            % In which a cluster is not standard normally distributed.
            % Take care of outliers first, those can be overlapping spikes,
            % that should not influence the decision.
            
            nD = self.P.nSubDimensions;
            maxID = max(cls);
            for i=1:nClasses
                candidateProjections = [];
                myId = cls(i);
                idx = classes==myId;
                [F pcs] = util.dimReductionPCA(pwSpikesX(idx,:), nD);
                for d = 1:nD
                    % Since clusters can be expected to be approximately standard 
                    % nomally distributed, we can test that on every single
                    % dimension individually. 
%                     [F bTooManyRemoved] = removeOutliers(fetX(:,d));
%                     if ~bTooManyRemoved && ~isSingleCluster(F);
                    if ~sorters.EMClusteringHirarchical.isSingleCluster(F(:,d));
                        candidateProjections = [candidateProjections d];
                    end
                end
                % Check if any cluster was thought to be not a nice single
                % cluster
                if isempty(candidateProjections)
                    continue
                end                
                myIDs  = clustering.gmmClu(F(:,candidateProjections),'kmin', 1,...
                                'kmax', 5,....
                                'EXPSIG', self.P.EXPSIG);
                myNewClasses = unique(myIDs);
                newIDs = (1:length(myNewClasses))+maxID;
                maxID = max(newIDs);
                for newClassIdx = 1:length(myNewClasses)
                    myIDs(myIDs==myNewClasses(newClassIdx)) = newIDs(newClassIdx);
                end
                classes(idx) = myIDs;
            end
        end
    end
    
    methods (Static)
        %%% ------------------------------------------------------
        function b = isSingleCluster(F)
            F = F - median(F);
            b = sum(abs(F)>10) < .10*length(F);
        end
        %%% ------------------------------------------------------
%         function [F bTooManyOutliers] = removeOutlier(F)
%             % Since F is supposed to be subspace of the prewhitened space
%             % Clusters should be standard normally distributed
%             if sum(abs(F)>10)>length(F)*.05
%                 F = [];
%                 bTooManyOutliers = true;
%                 return
%             else
%                 F = 
%             end
%         end
    end
end
    
    