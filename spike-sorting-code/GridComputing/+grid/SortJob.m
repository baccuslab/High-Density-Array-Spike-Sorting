%Example:
% sj = grid.SortJob(jobName, dataFiles, destFolder);
% sj.createBOTMGroups();
% sj.copyDataFiles();
% sj.setTaskParameters();
% sj.prepareTasks();
%
% --> Submit to grid
%
% all_tasks_completed = sj.waitForTasksToFinish();
% if all_tasks_completed
%     sj.copyBackResults();
%     sj.postprocessBOTMSorting();
%     sj.runQC()
% end

classdef SortJob < grid.GridJob
    properties (SetAccess=private)
    end
    
    properties
        sortJobP
        groupsidx
        MES
    end
    
    methods
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = SortJob(jobName, dataFiles, destFolder, varargin)
            self = self@grid.GridJob(jobName, varargin{:});
            self.taskType = 'SortJob';
            
            p = struct;
            p.useFilter = 0;
            p.nGroups = [];
            p.postProcOnly = false;
            
            p = mysort.util.parseInputs(p, self.P_untreated, 'error');
            self.sortJobP = p;
            
            %% Set the data files:
            if ~iscell(dataFiles) dataFiles = {dataFiles}; end
            self.startlocation.files.data = dataFiles;
            self.destinationlocation.folders.main = destFolder;
            
            %% Set the 'completed' variables to false:
            self.completed.createBOTMGroups = false;
            self.completed.setTaskParameters = false;
            self.completed.copyBackResults = false;
            self.completed.copyDataFiles = false;
            self.completed.postprocessBOTMSorting = false;
            
            %% Set the necessary parameters when only the postprocessing step needs to be done:
            if self.sortJobP.postProcOnly
                self.completed.createBOTMGroups = true;
                self.completed.setTaskParameters = true;
                self.completed.prepareTasks = true;
                self.completed.waitForTasksToFinish = true;
                self.completed.copyBackResults = true;
                self.completed.copyDataFiles = true;
                
                self.destinationlocation.files.groupFile = fullfile(self.scratchlocation.folders.main, 'groupFile.mat');
                
                %if ~isempty(self.P.dataPath)
                %    self.scratchlocation.folders.groups = fullfile(self.scratchlocation.folders.main);
                %else
                self.scratchlocation.folders.groups = fullfile(self.scratchlocation.folders.main, 'output');
                %end
                
                self.destinationlocation.files.groupPaths = self.findSubFolders(self.scratchlocation.folders.groups);
            end
        end
        
        % -----------------------------------------------------------------
        function createBOTMGroups(self)
            if self.completed.createBOTMGroups return; end
            
            disp('Create DataStructure object...');
            DSFull = mysort.mea.CMOSMEA(self.startlocation.files.data, 'useFilter', self.sortJobP.useFilter, 'name', self.jobName);%'PREFILT');
            self.MES = DSFull.MultiElectrode.toStruct();
            
            %% Create groups based on the electrode positions and save them to a file:
            self.scratchlocation.files.groupFile = fullfile(self.scratchlocation.folders.main, 'groupFile.mat');
            %if ~isempty(self.P.dataPath)
            %    self.destinationlocation.files.groupFile = self.scratchlocation.files.groupFile;
            %end
            
            try
                load(self.scratchlocation.files.groupFile)
                disp('Groups already exist')
            catch
                disp('Create groups...');
                electrodePositions = DSFull.MultiElectrode.electrodePositions;
                electrodeNumbers   = DSFull.MultiElectrode.electrodeNumbers;
                [groupsidx nGroupsPerElectrode] = mysort.mea.constructLocalElectrodeGroups(electrodePositions(:,1), electrodePositions(:,2), 'maxElPerGroup', 9);
                
                %% Limit number of groups if necessary:
                if ~isempty(self.sortJobP.nGroups)
                    groupsidx = {groupsidx{1:self.sortJobP.nGroups}};
                end
                
                disp(['Number of groups created: ' num2str(length(groupsidx))])
                
                groups = {};
                for ii= 1:length(groupsidx)
                    groups{ii} = electrodeNumbers(groupsidx{ii});
                end
                
                save(self.scratchlocation.files.groupFile, 'groups', 'electrodeNumbers', 'electrodePositions', 'nGroupsPerElectrode', 'groupsidx');
            end
            
            self.groupsidx = groupsidx;
            self.startIndex = 1;
            self.endIndex = length(self.groupsidx);
            self.taskIDs = self.startIndex:self.endIndex;
            
            
                %% Create SortJob specific output folder (on scratch):
                self.scratchlocation.folders.groups = fullfile(self.scratchlocation.folders.main, 'groups');
                [dir_exists,mess,messid] = mkdir( self.scratchlocation.folders.groups );
                assert(dir_exists, 'Output directory could not be created!');
            
            
            self.completed.createBOTMGroups = true;
        end
        
        function setTaskParameters(self)
            if self.completed.setTaskParameters return; end
            taskParameters = struct;
            
            %% Sorting parameters form ana.startHDSorting:
            taskParameters.sortingParameters = struct;
            taskParameters.sortingParameters.spikeDetection.method = '-';
            taskParameters.sortingParameters.spikeDetection.thr = 4.2;
            taskParameters.sortingParameters.artefactDetection.use = 0;
            taskParameters.sortingParameters.botm.run = 0;
            taskParameters.sortingParameters.spikeCutting.maxSpikes = 200000000000; % Set this to basically inf
            
            taskParameters.sortingParameters.noiseEstimation.minDistFromSpikes = 80;
            
            taskParameters.sortingParameters.spikeAlignment.initAlignment = '-';
            taskParameters.sortingParameters.spikeAlignment.maxSpikes = 50000;     % so many spikes will be clustered
            taskParameters.sortingParameters.clustering.maxSpikes = taskParameters.sortingParameters.spikeAlignment.maxSpikes;  % dont align spikes you dont cluster...
            taskParameters.sortingParameters.clustering.meanShiftBandWidth = sqrt(1.8*6);
            taskParameters.sortingParameters.mergeTemplates.merge = 1;
            taskParameters.sortingParameters.mergeTemplates.upsampleFactor = 3;
            taskParameters.sortingParameters.mergeTemplates.atCorrelation = .93; % DONT SET THIS TOO LOW! USE OTHER ELECTRODES ON FULL FOOTPRINT TO MERGE
            taskParameters.sortingParameters.mergeTemplates.ifMaxRelDistSmallerPercent = 30;
            %% ---------------------------------------------
            
            taskParameters.dataFiles = self.scratchlocation.files.data;
            taskParameters.runName = self.jobName;
            taskParameters.outputPath = self.scratchlocation.folders.groups;
            
            disp('Create task files from groups...');
            %% Create cell variable allTaskParamters:
            for ii = 1:length(self.taskIDs)
                
                self.scratchlocation.files.groupPaths{ii} = ['group' sprintf('%04d', self.taskIDs(ii))];
                
                taskParameters.groupidx = self.groupsidx{ii};
                taskParameters.taskID = self.taskIDs(ii);
                taskParameters.groupPath = self.scratchlocation.files.groupPaths{ii};
                
                self.allTaskParameters{ii} = taskParameters;
            end
            
            self.completed.setTaskParameters = true;
        end
        
        % -----------------------------------------------------------------
        function postprocessBOTMSorting(self)
            % Process all local sortings into a final sorting.
            if self.completed.postprocessBOTMSorting return; end
            
            disp('Start postprocessBOTMSorting...')
            self.destinationlocation.files.results = fullfile(self.destinationlocation.folders.main, [self.jobName '_results.mat']);
            try
                load(self.destinationlocation.files.results)
                disp('Postprocessing already finished.')
            catch
                disp('Load group information...');
                GF = load(self.destinationlocation.files.groupFile, 'groups', 'electrodeNumbers', 'electrodePositions', 'nGroupsPerElectrode', 'groupsidx');
                
                disp('Start postprocessing...');
                [gdf_merged, T_merged, localSorting, localSortingID, G] = mysort.HDSorting.processLocalSortings(...
                    self.destinationlocation.folders.groups,...
                    self.jobName, GF.groups, GF.groupsidx, ...
                    self.destinationlocation.files.groupPaths);
                
                disp('Check and save data...');
                units = unique(gdf_merged(:,1));
                nU = length(units);
                assert(length(localSorting) == nU, 'must be identical');
                assert(length(localSortingID) == nU, 'must be identical');
                assert(size(T_merged,3) == nU, 'must be identical');
                
                disp('Saving postprocessing results...')
                save(self.destinationlocation.files.results, 'gdf_merged', 'T_merged', 'localSorting', 'localSortingID', 'G', '-v7.3');
            end
            
            self.completed.postprocessBOTMSorting = true;
        end
        
        
        %         % -----------------------------------------------------------------
        %         function createSummaryFile(self)
        %
        %             self.destinationlocation.files.summary = fullfile(self.destinationlocation.folders.main, 'summary.mat');
        %             try
        %                 load(self.destinationlocation.files.summary)
        %             catch
        %                 try
        %                     res = load(self.destinationlocation.files.results)
        %
        %                     units = unique( res.gdf_merged(:,1) );
        %                     nUnits = length(units);
        %                     %[dir_exists,mess,messid] = mkdir(self.destinationlocation.folders.main, 'qcplots');
        %                     %self.destinationlocation.folders.qcplots = fullfile( self.destinationlocation.folders.main, 'qcplots');
        %
        %                     %DSFull = mysort.mea.CMOSMEA(self.startlocation.files.data);%, 'useFilter', self.sortJobP.useFilter, 'name', self.jobName);
        %                     %MES = DSFull.MultiElectrode.toStruct();
        %                     save(self.destinationlocation.files.summary, 'units', 'nUnits');
        %
        %                 catch
        %                     disp('Creation of summary file failed!');
        %                     return;
        %                 end
        %             end
        %         end
        
        
        % -----------------------------------------------------------------
        function runQC(self)
            try
                res = load(self.destinationlocation.files.results)
                
                [dir_exists,mess,messid] = mkdir(self.destinationlocation.folders.main, 'qcplots');
                self.destinationlocation.folders.qcplots = fullfile( self.destinationlocation.folders.main, 'qcplots');
                
                DSFull = mysort.mea.CMOSMEA(self.startlocation.files.data);%, 'useFilter', self.sortJobP.useFilter, 'name', self.jobName);
                MES = DSFull.MultiElectrode.toStruct();
                
            catch
                disp('Quality control failed!');
                return;
            end
            
            figures = struct;
            
            %% ISIH:
            self.destinationlocation.files.isih = fullfile( self.destinationlocation.folders.qcplots, 'isih');
            if ~exist(self.destinationlocation.files.isih, 'file')
                F = mysort.plot.isi( res.gdf_merged )
                mysort.plot.savefig(F.figureHandle, self.destinationlocation.files.isih)
            end
            
            %% Footprints whole
            self.destinationlocation.files.footprints = fullfile( self.destinationlocation.folders.qcplots, 'footprints_whole');
            if ~exist(self.destinationlocation.files.footprints, 'file')
                F.figureHandle = figure();
                mysort.plot.waveforms2D(res.T_merged, MES.electrodePositions, 'IDs', res.localSortingID);
                mysort.plot.savefig(F.figureHandle, self.destinationlocation.files.footprints)
            end
            
            %% Footprints localized
            self.destinationlocation.files.footprints2 = fullfile( self.destinationlocation.folders.qcplots, 'footprints_localized');
            if ~exist(self.destinationlocation.files.footprints2, 'file')
                F.figureHandle = figure(); P.AxesHandle = []
                for i = 1:length(res.localSorting)
                    id = res.localSortingID(i);
                    lu = res.localSorting(i);
                    P = mysort.plot.waveforms2D(0.1*res.T_merged(:,:,i), MES.electrodePositions, 'IDs', (1000*id+lu), 'maxNumberOfChannels', 10, 'AxesHandle', P.AxesHandle, 'plotArgs', {'color', mysort.plot.vectorColor(i)});
                end
                mysort.plot.savefig(F.figureHandle, self.destinationlocation.files.footprints2)
            end
            
        end
        
        % -----------------------------------------------------------------
        function copyBackResults(self, copyEverything)
            % Copy results from the grid folder back to the destination folder
            if self.completed.copyBackResults return; end
            
            if ~isempty(self.P.dataPath) 
                self.destinationlocation.files.groupFile = self.scratchlocation.files.groupFile;
                self.destinationlocation.files.summary = self.scratchlocation.files.summary;
                self.destinationlocation.folders.groups = self.scratchlocation.folders.groups;
                self.destinationlocation.files.groupPaths = self.findSubFolders(self.scratchlocation.folders.groups);
                self.destinationlocation.folders.main = self.scratchlocation.folders.main; % sorting output here
                self.destinationlocation.folders.groups = self.scratchlocation.folders.groups; % sorting output here
                return;
            end
            
            if nargin < 2 copyEverything = false; end
            
            disp('Copy back result file(s)...');
            
            [dir_exists,mess,messid] = mkdir(self.destinationlocation.folders.main); %self.folders.destination);
            assert(dir_exists, 'Destination directory could not be created!');
            
            self.destinationlocation.files.groupFile = self.copyFile(self.scratchlocation.files.groupFile, self.destinationlocation.folders.main);
            self.destinationlocation.files.summary = self.copyFile(self.scratchlocation.files.summary, self.destinationlocation.folders.main);
            
            self.destinationlocation.folders.groups = fullfile( self.destinationlocation.folders.main, 'groups');
            [dir_exists,mess,messid] = mkdir(self.destinationlocation.folders.groups); %self.folders.destination);
            assert(dir_exists, 'Destination groups directory could not be created!');
            
            if copyEverything
                disp('Copy back entire sorting output...')
                self.destinationlocation.folders.groups = self.copyFolderContent(self.scratchlocation.folders.groups, self.destinationlocation.folders.groups);
                
                disp('Copy back data files...')
                self.destinationlocation.folders.data = self.copyFolder(self.scratchlocation.folders.data, self.destinationlocation.folders.main);
                disp('Done copying.')
            else
                %% Only copy back necessary (and small) output files:
                filesToCopy = {'_templates.mat', ...
                    '.P.mat',...
                    '.030spikes_det_merged.mat', ...
                    '.060cov.mat', ...
                    '.090clusters_meanshift.mat',...
                    '.100botm_matching.mat', ...
                    '.110clusters_meanshift_merged.mat'};
                
                groupFolders = self.findSubFolders(self.scratchlocation.folders.groups)
                
                for i = 1:length(groupFolders)
                    groupFolder = groupFolders{i};
                    [pathstr,name,ext] = fileparts(groupFolder);
                    newGroupFolder = fullfile(self.destinationlocation.folders.groups, name);
                    
                    [dir_exists,mess,messid] = mkdir(newGroupFolder);
                    assert(dir_exists, 'Destination group folder could not be created!');
                    
                    for j = 1:length(filesToCopy)
                        fullFile = fullfile(self.scratchlocation.folders.groups, groupFolder, [self.jobName filesToCopy{j}])
                        self.copyFile(fullFile, newGroupFolder);
                    end
                    
                    self.destinationlocation.files.groupPaths{i} = groupFolder; %newGroupFolder;
                end
            end
            
            disp('Sorting results have been copied to:')
            disp(self.destinationlocation.folders.main)
            
            self.completed.copyBackResults = true;
        end
        
        
        % -----------------------------------------------------------------
        function copyDataFilesToBinary(self)
            % This function copies the original h5 file to a binary file and a
            % h5 metadata file in the job folder. From now on, meta-file will
            % be treated as the original data file.
            % The binary files are created directly on the scratch.
            disp('Copy data file(s) to binary...');
            for i = 1:length(self.startlocation.files.data)
                
                [pathstr,name,ext] = fileparts(self.startlocation.files.data{i});
                
                fileNameH5 = fullfile(self.scratchlocation.folders.data, [name '_meta' ext]);
                fileNameBin = fullfile(self.scratchlocation.folders.data, [name '.dat']);
                
                
                if exist(fileNameBin, 'file') ~= 2
                    mysort.mea.copyH5toBinary(self.startlocation.files.data{i}, fileNameH5, fileNameBin);
                end
                
                %% Remember the location of the datafiles on the scratch now:
                self.scratchlocation.files.data{i} = fileNameH5;
                
                disp([ num2str(i) ' of ' num2str(length(self.startlocation.files.data)) ' copied.']);
            end
        end
        
        function copyDataFiles(self)
            if self.completed.copyDataFiles return; end
            
            if ~isempty(self.P.dataPath)
                self.scratchlocation.files.data = self.startlocation.files.data;
                return;
            end
            
            disp('Copy data file(s)...');
            for i = 1:length(self.startlocation.files.data)
                
                %% Check if the data (the first file at least) is binary:
                m_test = mysort.mea.CMOSMEA(self.startlocation.files.data{i});
                
                if ~m_test.isBinaryFile()
                    [pathstr,name,ext] = fileparts(self.startlocation.files.data{i});
                    fileNameH5 = fullfile(self.scratchlocation.folders.data, [name '.h5']);
                    fileNameBin= fullfile(self.scratchlocation.folders.data, [name '.dat']);
                    
                    %% Create a binary file on the scratch:
                    if exist(fileNameBin, 'file') ~= 2
                        mysort.mea.copyH5toBinary(self.startlocation.files.data{i}, fileNameH5, fileNameBin);
                    end
                    self.scratchlocation.files.data{i} = fileNameH5;
                else
                    [pathstr,name,ext] = fileparts(self.startlocation.files.data{i});
                    fileNameH5 = self.startlocation.files.data{i}; %fullfile(self.startlocation.folders.data, [name '_meta' ext]);
                    fileNameBin = fullfile(pathstr, [name '.dat']);
                    
                    %% Copy meta and binary file:
                    self.scratchlocation.files.data{i} = self.copyFile(fileNameH5, self.scratchlocation.folders.data) ;
                    self.copyFile(fileNameBin,self.scratchlocation.folders.data);
                end
                
                disp([ num2str(i) ' of ' num2str(length(self.startlocation.files.data)) ' copied.']);
            end
            self.completed.copyDataFiles = true;
        end
        
        %% Auxilary function for finding non-hidden subfolders
        function list = findSubFolders(self, folder)
            d = dir(folder);
            isub = [d(:).isdir]; %# returns logical vector
            list = {d(isub).name}';
            list(ismember(list,{'.','..'})) = [];
        end
    end
    
    methods(Static)
        %------------------------------------------------------------------
        %------------------------- RUN  FUNCTION ---------------------------
        %------------------------------------------------------------------
        function run(taskFile, debugFlag)
            if nargin < 2
                debugFlag = false;
            end
            
            %% List of default parameters:
            sortP = struct;
            sortP.botm.Tf = 55;
            sortP.botm.cutLeft = 10;
            sortP.botm.run = 0;
            sortP.spikeCutting.blockwise = false;
            
            taskP = struct;
            taskP.runName = 'BOTMdefault';
            
            %% Load taskFile:
            T = load(taskFile);
            
            %% Write "taskParameters" to struct "taskP" and "sortP":
            sortP = mysort.util.mergeStructs(sortP, T.taskParameters.sortingParameters);
            taskParameters = rmfield(T.taskParameters,'sortingParameters');
            
            taskP = mysort.util.mergeStructs(taskP, T.taskParameters);
            clear T;
            
            %% Check necessary parameters:
            assert( isfield(taskP, 'dataFiles'), 'Task aborted: field taskParameters.dataFiles not specified!');
            for f = 1:length(taskP.dataFiles)
                % Note: exist(...) == 2: name is a full pathname to a file
                assert( exist(taskP.dataFiles{f}, 'file') == 2, ['Task aborted: dataFile ' taskP.dataFiles{f} ' not found!']);
            end
            assert( exist(taskP.outputPath, 'dir') == 7, ['Task aborted: no valid outputPath specified! Path given:' taskP.outputPath]);
            assert( isfield(taskP, 'groupidx'), 'Task aborted: field taskParameters.groupidx not specified!');
            assert( isfield(taskP, 'taskID'), 'Task aborted: field taskParameters.taskID not specified!');
            
            %% (Re-)Set reporting file:
            rep = mysort.ds.binaryFileMatrix(taskP.reportFile, [1 2], 'writable', true);
            rep(:,:) = [0 0];
            
            if ~debugFlag
                try
                    mainBlock();
                catch ME
                    errorHandling(ME);
                end
            else
                mainBlock();
            end
            
            function mainBlock()
                disp('Creating DS...')
                %% Sort:
                DS = mysort.mea.CMOSMEA(taskP.dataFiles, 'useFilter', 0, 'name', 'PREFILT');
                MES = DS.MultiElectrode.toStruct();
                
                DS.restrictToChannels(taskP.groupidx);
                disp('Start sorting...')
                [S P_] = mysort.sorters.sort(DS, fullfile(taskP.outputPath, taskP.groupPath ), taskP.runName, sortP);
                
                %% RELEASE CHANNEL RESTRICTIONS FOR TEMPLATE ESTIMATION
                DS.restrictToChannels();
                disp('Starting template estimation...')
                mysort.HDSorting.startHDSortingTemplateEstimation(taskP.outputPath, taskP.groupPath, taskP.runName, sortP.botm.Tf, sortP.botm.cutLeft, DS, taskP.groupidx, MES);
                
                %% Write to reporter file:
                disp('Writing results...')
                rep = mysort.ds.binaryFileMatrix(taskP.reportFile, [1 2], 'writable', true);
                rep(:,:) = [1 0];
            end
            
            function errorHandling(ME)
                
                disp('Catch error...')
                errStr = mysort.util.buildLastErrString();
                disp(errStr)
                
                rep = mysort.ds.binaryFileMatrix(taskP.reportFile, [1 2], 'writable', true);
                rep(:,:) = [0 1];
                rethrow(ME)
            end
        end
        
    end
end
