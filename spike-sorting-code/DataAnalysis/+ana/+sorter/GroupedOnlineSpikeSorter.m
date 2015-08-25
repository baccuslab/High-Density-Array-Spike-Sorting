classdef GroupedOnlineSpikeSorter < handle
    properties (SetAccess=private)
        
    end
    properties
        P
        SD
        group_gdfs
        gdfs
        singleElGroupSorters
        electrodeGroups
        savePath
        buffer
        
        trialNames 
        trialCounter
        chunkCounter
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = GroupedOnlineSpikeSorter(electrodeGroups, savePath, varargin)
            self.P.mergeSpikesMaxDist = 20;
            self.P = mysort.util.parseInputs(self.P, varargin, 'error');
            
            self.savePath = savePath;
            self.electrodeGroups = electrodeGroups;
            
            self.SD = ana.sorter.IndividualChannelSpikeDetector(...
                      'thresholdFactor', 4.5, ...
                      'minDist', 20);
            self.singleElGroupSorters = ana.sorter.OnlineSpikeSorter.empty(0, length(electrodeGroups));
            for i=1:length(electrodeGroups)
                group_savePath = fullfile(savePath, sprintf('group%03d', i));
                self.singleElGroupSorters(i) = ana.sorter.OnlineSpikeSorter(group_savePath);
            end
            
            self.buffer = {};
            self.chunkCounter = 0;
            self.trialCounter = 0;
            self.trialNames = {};
            self.group_gdfs = cell(0, length(electrodeGroups));
            self.gdfs = cell(0,1);
        end
        % --------------------------------------------------------
        function startNextTrial(self, trialName)
            self.trialNames{end+1} = trialName;
            self.trialCounter = self.trialCounter+1;
        end
        
        % --------------------------------------------------------
        function processChunk(self, X)
            % create handle class to avoid data copying
%             X = mysort.ds.Matrix(X);            
            
            self.SD.processChunk(X);
            
            if self.SD.bIsTrained
                if ~isempty(self.buffer)
                    for i=1:length(self.buffer)
                        self.processChunk_(self.buffer{i});
                        self.buffer{i} = [];
                    end
                    self.buffer = {};
                end
                
                self.processChunk_(X);
            else
                self.buffer{end+1} = X;
                return
            end
        end
        % --------------------------------------------------------
        function processChunk_(self, X)
            self.chunkCounter = self.chunkCounter+1;
            
            nG = length(self.electrodeGroups);
            G = cell(1, nG);
            for g = 1:nG
                group = self.electrodeGroups{g};
                S = struct();
                S.detUp = {};
                S.pksUp = {};
                S.detDown = stdown(group)';
                S.pksDown = pksdown(group)';
                S.X = X(:, group);             
                G{g} = S;
            end   
            
            spso = self.singleElGroupSorters;
            for g=1:nG
                allspikes = ana.sorter.mergeSpikes(G{g}.detUp, G{g}.pksUp, G{g}.detDown, G{g}.pksDown, self.P.mergeSpikesMaxDist);
                spso(g).processChunk(G{g}.X, allspikes(:,1), chunkName);
            end
        end
    end
end
% 
% 
% %     for i=1:length(flist)
% %         myFile = fullfile(dpath, flist{i});
% %     for g=1:nG
% %         gdfg = [[CLU(:).ids]
% %     
% %             
% %         gdf = 
% %     end
% 
% 
%     if 0
%         G = ana.mergeLocalSortings(G, meanNoiseStd); 
% %         cov(CLU{g}.F(CLU{g}.ids==6,:))        
%     end
%     if 0
% 
%         g = 18;
%         figure;
%         ah = mysort.plot.subplot2([1 3]);
%         mysort.plot.clustering(CLU{g}.F', CLU{g}.ids, [], [], 'axesHandles', ah(2:3));
%         mysort.plot.waveformsVertical(CLU{g}.templates, 'IDs', CLU{g}.classes, ...
%             'linewidth', 2, 'channelSpaceFactor',1, 'axesHandle', ah(1))
%         
%         mysort.plot.clustering(CLU{g}.F', CLU{g}.ids_clu, [], []);
%     end
%     
%     if 0 
%         %%
%         g=1;
%         vT = mysort.wf.t2v(CLU{g}.templates);
%         nT = size(vT,1);
%         CC = zeros(nT,nT);
%         for i=1:nT
%             for j=i+1:nT
%                 cc = corrcoef(vT([i j],:)');
%                 CC(i,j) = cc(2,1);
%             end
%         end
%         figure;
%         subplot(1,2,1)
%         imagesc(CC)
%         colorbar
%         subplot(1,2,2)
%         imagesc(CC>.9)
%         colorbar        
%     end
%     
%     