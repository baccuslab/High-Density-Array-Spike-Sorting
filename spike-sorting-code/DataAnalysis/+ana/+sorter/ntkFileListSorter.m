classdef ntkFileListSorter < mysort.util.HandleObject
    properties (SetAccess=private)
        
    end
    properties
        P
        dataPath
        savePath
        fileList
        fileDataBuffer
        
        electrodes
        ntkW_first
        SD
        singleElGroupSorters 
        gdfs
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = ntkFileListSorter(dataPath, fileList, savePath, varargin)
            self.P.minSamplesForTraining = 10*60*20000; % in min
            self.P.mergeSpikesMaxDist = 20;
            self.P = mysort.util.parseInputs(self.P, varargin, 'error');
            
            self.dataPath = dataPath;
            self.fileList = fileList;
            self.savePath = savePath;
            
            self.init();

%             self.templateMatching();
        end
        % --------------------------------------------------------
        function save(self)
            saveFile = fullfile(self.savePath, 'ntkFileSorterSave.mat');
            S = self.toStruct();
            for i=1:length(self.singleElGroupSorters)
                se(i) = self.singleElGroupSorters(i).toStruct();
            end
            S.singleElGroupSorters = se;
            clear se
            save(saveFile, 'S');
        end
        % --------------------------------------------------------
        function V = checkFList(self)
            fList = self.fileList;
            V = zeros(2, length(fList));
            % Walk over .ntk files, check for existance and ME
            for i=1:length(fList)
                % init file and load data
                myFile = fullfile(self.dataPath, fList{i});
                if ~exist(myFile, 'file')
                    disp(['Not Found: ' fList{i}])
                    V(:, i) = [-1 -1]';
                    continue
                end
                s = dir(myFile);
                V(1,i) = s.bytes;
                ntkW = ana.sorter.ntkFileWrapper(myFile);
                if self.ntkW_first.MultiElectrode == ntkW.MultiElectrode
                    V(2,i) = 1;
                end
            end
            if any(V(1,:)==-1)
                disp('There were file not found errors!');
            end
            if any(V(2,:)<1)
                disp('There were non matching configurations!');
            end
            if all(V(:)>0)
                disp('##### Flist is fine!! #####');
            end
        end        
        % --------------------------------------------------------
        function load(self)
            saveFile = fullfile(self.savePath, 'ntkFileSorterSave.mat');
            L = load(saveFile);
            self.fromStruct(L.S);
            for i=1:length(self.singleElGroupSorters)
                se(i) = ana.sorter.OfflineSpikeSorter();
                se(i).fromStruct(self.singleElGroupSorters(i));
            end
            self.singleElGroupSorters = se;
        end        
        % --------------------------------------------------------
        function init(self)
            % Make groupings of electrodes to be sorted independently
            % use first file in flist to initialize sorting groups
            firstFile = fullfile(self.dataPath, self.fileList{1});
            self.ntkW_first = ana.sorter.ntkFileWrapper(firstFile);
            [groupsidx nGroupsPerElectrode] = self.ntkW_first.MultiElectrode.getLocalElectrodeGroups();
            self.SD = ana.sorter.IndividualChannelSpikeDetector(...
                      'thresholdFactor', 4.5, ...
                      'minDist', 20);
            self.electrodes.ME = self.ntkW_first.MultiElectrode.toStruct();
            self.electrodes.groups = groupsidx;
            self.electrodes.nGroupsPerElectrode = nGroupsPerElectrode;
            
            self.singleElGroupSorters = ana.sorter.OfflineSpikeSorter.empty(0, length(self.electrodes.groups));
            
            for i=1:length(self.electrodes.groups)
                self.singleElGroupSorters(i) = ana.sorter.OfflineSpikeSorter();
            end            
        end
        % --------------------------------------------------------
        function loadDataForTraining(self)
            fList = self.fileList;
            self.fileDataBuffer = cell(0,1);
            cumLength = 0;
            % Walk over .ntk files until trained
            for i=1:length(fList)
                chunkBuffer = cell(0,1);
                % init file and load data
                myFile = fullfile(self.dataPath, fList{i});
                ntkW = ana.sorter.ntkFileWrapper(myFile);
                assert(self.ntkW_first.MultiElectrode == ntkW.MultiElectrode, 'This ntk file has a different configuration!');
                while ~ntkW.eof()
                    chunkBuffer{end+1,1} = ntkW.getNextChunkFilteredData();
                    cumLength = cumLength + size(chunkBuffer{end,1},1);
                    if cumLength > self.P.minSamplesForTraining
                        self.fileDataBuffer{i,1} = cell2mat(chunkBuffer);
                        return
                    end
                end
                self.fileDataBuffer{i,1} = cell2mat(chunkBuffer);
            end            
        end
      
        % --------------------------------------------------------
        function train(self)
            for i=1:size(self.fileDataBuffer,1)
                self.SD.processChunk(self.fileDataBuffer{i,1});
            end
            
            % Get spikes
            nG = length(self.electrodes.groups);            
            cl = self.singleElGroupSorters(1).P.spikeCutting.cutLeft;
            cr = self.singleElGroupSorters(1).P.spikeCutting.Tf - cl +1;
            for i=1:size(self.fileDataBuffer,1)
                [stdown pksdown stup pksup] = self.SD.getSpikeTimesForChunk(i);
                
                % Process groups
                for g = 1:nG
                    group = self.electrodes.groups{g};
                    X = self.fileDataBuffer{i,1}(:, group); 

                    allspikes = ana.sorter.mergeSpikes({}, {}, stdown(group)', pksdown(group)', self.P.mergeSpikesMaxDist);
                    % remove spikes at the borders of chunk
                    allspikes(allspikes(:,1) <= cl,:) = [];
                    allspikes(allspikes(:,1) >= size(X,1)-cr,:) = [];                    
                    self.singleElGroupSorters(g).processChunk(X, allspikes(:,1));
                end                    
            end
            self.forceTrain();
            self.fileDataBuffer = {};
        end  
        
        % --------------------------------------------------------
        function forceTrain(self)
            % Check that all gorups are trained
            for i=1:length(self.singleElGroupSorters)
                if ~self.singleElGroupSorters(i).bIsTrained
                    fprintf('Group %d had not enough spikes, forcing to train!\n', i)
                    self.singleElGroupSorters(i).train();
                end
            end        
        end
        % --------------------------------------------------------
        function templateMatching(self)
            % check if trained
            for g=1:length(self.singleElGroupSorters)
                if ~self.singleElGroupSorters(g).bIsTrained
                    fprintf('A single El Sorter (%d) is not trained yet. Aborting. Call train first!\n', g);
                    return
                end
            end
            gdfs = cell(0,1);
            startTrial = 1;            
            if ~isempty(self.gdfs)
                gdfs = self.gdfs;
                startTrial = size(gdfs,1)-1;
            end
            fList = self.fileList;
            nG = length(self.electrodes.groups);    
            cl = self.singleElGroupSorters(1).P.spikeCutting.cutLeft;
            cr = self.singleElGroupSorters(1).P.spikeCutting.Tf - cl +1;
            % Walk over .ntk files until trained
            for i=startTrial:length(fList)
                chunkBuffer = cell(0,nG);
                chunk = 1;
                chunkTimeOffset = 0;
                % init file and load data
                myFile = fullfile(self.dataPath, fList{i});
                ntkW = ana.sorter.ntkFileWrapper(myFile);
                assert(self.ntkW_first.MultiElectrode == ntkW.MultiElectrode, 'This ntk file has a different configuration!');
                while ~ntkW.eof()
                    X = ntkW.getNextChunkFilteredData();
                    if size(X,1) < 200
                        break;
                    end
                    self.SD.processChunk(X);
                    [stdown pksdown stup pksup] = self.SD.getSpikeTimesForChunk(self.SD.chunkCounter);
                    % Process groups
                    for g = 1:nG
                        group = self.electrodes.groups{g};

                        allspikes = ana.sorter.mergeSpikes({}, {}, stdown(group)', pksdown(group)', self.P.mergeSpikesMaxDist);
                        % remove spikes at the borders of chunk
                        allspikes(allspikes(:,1) <= cl,:) = [];
                        allspikes(allspikes(:,1) >= size(X,1)-cr,:) = [];                    
                        
                        if isempty(allspikes)
                            chunkBuffer{chunk, g} = [];
                        else
                            self.singleElGroupSorters(g).processChunk(X(:, group), allspikes(:,1));
                            chunkGroupGdf = self.singleElGroupSorters(g).gdfs{end};
                            chunkGroupGdf(:,2) = chunkGroupGdf(:,2) + chunkTimeOffset;
                            chunkBuffer{chunk, g} = chunkGroupGdf;                            
                        end

                    end    
                    chunkTimeOffset = chunkTimeOffset + size(X,1);
                    chunk = chunk +1;                    
                end
                mygdfs = cell(1,nG);
                for g=1:nG
                    mygdfs{1,g} = cell2mat(chunkBuffer(:,g));
                end
                self.gdfs{i,1} = mygdfs;
            end    
%             self.gdfs = gdfs;
        end          
    end
end