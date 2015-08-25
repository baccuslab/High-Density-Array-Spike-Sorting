%% INPUT CHECKS
assert(exist('CONFIG_FILE_LOCATION', 'var')>0, 'Must be set');
assert(exist('submitForReal', 'var')>0, 'Must be set');
assert(exist('runLocally', 'var')>0, 'Must be set');
assert(exist('bCorrectForMissingFrames', 'var')>0, 'Must be set');
assert(exist('INTERMEDIATE_DATA_PATH', 'var')>0, 'Must be set');
assert(exist('expNames', 'var')>0, 'Must be set');

%% PREPROCESSING AND SUBMISSION TO GRID ENGINE
t_start_all = tic;
clear PATHDEFS;

SJ = cell(1, length(expNames));
all_Experiments_H5Files = {};
for eidx = 1:length(expNames)
    expName = expNames{eidx};
    fprintf('Preprocessing %s\n', expName);
    
    %% INIT AND CREATE FOLDER NAMES
    sortingName = 'Sorting1';
    
    pd = pdefs();
    clear pathdefs;
    % Make sure the jobname does not start with a number
    pathdefs.jobName =  ['J' expName];
    pathdefs.inputPath = CONFIG_FILE_LOCATION;
    pathdefs.preprocessedOutPath = fullfile(INTERMEDIATE_DATA_PATH, expName, 'Preprocessed');
    pathdefs.sortingOutPath      = fullfile(INTERMEDIATE_DATA_PATH, expName, 'SortOut');
    pathdefs.exportResultPath    = fullfile(INTERMEDIATE_DATA_PATH, expName, 'SortResult');
    
    
    %% PREPROCESS MEA1K FILES
    
    if ~exist(pathdefs.preprocessedOutPath, 'file'); mkdir(pathdefs.preprocessedOutPath); end
    confListFile = fullfile(pathdefs.inputPath, [expName '.mat']);
    c = load(confListFile);
    pathdefs.configNTKLists = c.configNTKLists;
    
    % HACK BECAUSE MICHELE MOVED THE ROSKA FOLDER ONE UP
    for i=1:length(pathdefs.configNTKLists)
        for j=1:length(pathdefs.configNTKLists{i})
            pathdefs.configNTKLists{i}{j} = strrep(pathdefs.configNTKLists{i}{j}, '/recordings/Mea1k/michelef/Roska/', '/recordings/Mea1k/michelef/');     
        end
    end
    % END HACK
            
    bIsMea1Krecording = true;
    try
        nrk_files_per_config = c.xy_pos;
        % HACK BECAUSE MICHELE MOVED THE ROSKA FOLDER ONE UP
        for j=1:length(nrk_files_per_config)
            if iscell(nrk_files_per_config{j})
                for jj=1:length(nrk_files_per_config{j})
                    nrk_files_per_config{j}{jj} = strrep(nrk_files_per_config{j}{jj}, '/recordings/Mea1k/michelef/Roska/', '/recordings/Mea1k/michelef/');     
                end
            else
                nrk_files_per_config{j} = strrep(nrk_files_per_config{j}, '/recordings/Mea1k/michelef/Roska/', '/recordings/Mea1k/michelef/');     
            end
        end     
        % END HACK
    catch
        % THIS IS AN OLD HIDENS RECORDING
        bIsMea1Krecording = false;
        nrk_files_per_config = [];
    end

    %% Now we CONVERT all ntk files to H5. Doing this we also prefilter the data
    allFiles = {};
    allCombis = [];
    pathdefs.h5FileLists = {};
    mapPerFile = {};
    for i=1:length(pathdefs.configNTKLists)
        if bIsMea1Krecording
            % load config file
            fname = nrk_files_per_config{i};
            if iscell(fname)
                fname = fname{1};
            end

            fname = strrep(fname, '.u.', '.');
            try 
                X = dlmread(fname);
            catch
                fname = strrep(fname, 'config/config_', 'config/configs_');
                try
                    X = dlmread(fname);
                catch
                    fname = strrep(fname, 'config/config', 'config/config_');
                    X = dlmread(fname);
                end
            end
            isNotConnected = any(X==-1,2);
            X(isNotConnected,:) = [];
            % Im H5 file ist spalte 0 channel 0.
            % 1023 ist der letzte Kanal dann DAC und !command counter"
            % DAC ist nicht im nrk file 
            idx = X(:,1) + 1;
            m = struct();
            m.mposx(idx) = X(:,2);
            m.mposy(idx) = X(:,3);
            m.elNo(idx)   = X(:,1);
            m.chans(idx)  = X(:,1);
            m.noChans = length(m.elNo);
            m.file = fname;
        end
        for k=1:length(pathdefs.configNTKLists{i})
            fn = pathdefs.configNTKLists{i}{k};
            allFiles{end+1} = fn;
            if bIsMea1Krecording
                mapPerFile{end+1} = m;
            end
            [ffolder ffname fsuffix] = fileparts(fn);
            fnpf = fullfile(pathdefs.preprocessedOutPath, [ffname fsuffix '_prefilter.h5']);
            pathdefs.h5FileLists{i}{k} = fnpf;
        end
        all_Experiments_H5Files{eidx} = pathdefs.h5FileLists;
    end

    tic
%     disp('Combi 0')
    nF = length(allFiles);
    parfor combi=1:nF
%         disp(nF)
%         try
        disp(combi)
%         disp('Combi 1')
        fn = allFiles{combi};
%         disp('Combi 2')
        [ffolder ffname fsuffix] = fileparts(fn);
        fnpf = fullfile(pathdefs.preprocessedOutPath, [ffname fsuffix '_prefilter.h5']);
        bConvertFile = true;
        disp(['Checking if file exists ' num2str(combi)])
        if exist(fnpf, 'file')
            try
                bIsInProcess = hdf5read(fnpf, '/bFileIsInProcess');
    %                 bIsInProcess = 1;
                if bIsInProcess
                    % This means the file is still in process. We assume no
                    % other user/matlab instance is writing so it is
                    % probably a left over from some other run and we
                    % delete it.
                    delete(fnpf); 
                    try
                        % In case this was a binary file, delete this as
                        % well
                        [a,b,c] = fileparts(fnpf);
                        delete([a b '.dat']);
                    catch
                    end
                else
                    % The file was already properly processed, we dont need
                    % to do anything.
                    bConvertFile = false;
                end
            catch
                delete(fnpf);
            end
        end
        if bConvertFile
            % Only convert file if it was non existent or deleted
            if bIsMea1Krecording
                disp(['Converting Mea1K ' num2str(combi)])
                mysort.mea.preprocessMea1kH5File(fn, mapPerFile{combi}, 'prefilter', 1, 'outFile', fnpf, 'h5path', h5Path);
            else
                disp(['Converting Highdens ' num2str(combi)])
                mysort.mea.convertNTK2HDF(fn, 'prefilter', 1, 'outFile', fnpf);
            end
        end
%         catch
%             errStr = mysort.util.buildLastErrString();
%             disp(errStr)
%         end
    end
    fn = [];
    toc
    
    disp('Done preprocessing!')
    
    %% START GRID PREPROCESSING
    pathdefs.sortOutPaths = {};
    SJ{eidx} = grid.SortJob.empty();
    nConfigs = length(pathdefs.h5FileLists);
    for i=1:nConfigs
        pathdefs.sortOutPaths{i} = fullfile(pathdefs.sortingOutPath, sortingName, ['Config' num2str(i)]);
        pathdefs.logFolders{i}   = fullfile(pathdefs.sortOutPaths{i}, 'log');
        sj = grid.SortJob(pathdefs.jobName, pathdefs.h5FileLists{i}, pathdefs.sortOutPaths{i}, ...
            'dataPath',  pathdefs.sortOutPaths{i},...
            'logFolder', pathdefs.logFolders{i});  % 'postProcOnly', postProcessingOnly);
        sj.createBOTMGroups();
        SJ{eidx}(i) = sj;
        SJ{eidx}(i).copyDataFiles();
        SJ{eidx}(i).setTaskParameters();
        SJ{eidx}(i).prepareTasks();  
        
        %% CREATE SUBMIT TOKENS ON SCRATCH TO TELL THE SUBMIT DEAMON TO SEND THE JOBS TO THE QUEUE
        if submitForReal % && i==2
            token_path = '/links/grid/scratch/frankef/tokens/';
            project_token = fullfile(token_path, ['start_' pathdefs.jobName '.tok']);
            sortOutPath = pathdefs.sortOutPaths{i};
            jobName = pathdefs.jobName;
            save(project_token, 'pathdefs', 'jobName', 'sortOutPath', 'expName');
            clear jobName;
            disp('Submitting Job to Grid-engine...');
            count = 0;
            bSubmitted = false;
            while ~bSubmitted && count < 120
                if mod(count-1, 5) == 0
                    fprintf('Waiting for job to be submitted... (%d sec)\n', count)
                end
                bSubmitted = ~exist(project_token, 'file');
                count = count+1;
                pause(1)
            end
            assert(bSubmitted, 'Token is still there, deamon dead?!');
            disp('Submit sucessful.');
        elseif runLocally
            sj_ = SJ{eidx}(i);
            for k=1:length(sj_.taskIDs)
                sj_.runTaskLocally(sj_.taskIDs(k));
            end
        end
    end
    if ~submitForReal
        disp('Not submitted to Grid, but all is prepared!');
    end
    PATHDEFS(eidx) = pathdefs;
end            

%% POSTPROCESSING AFTER GRID ENGINE IS DONE
parfor eidx = 1:length(expNames) %
    expName = expNames{eidx};  
    pathdefs = PATHDEFS(eidx);
    fprintf('Postprocessing %s\n', expName);
    nConfigs = length(SJ{eidx});
    
    %% WAIT FOR JOBS TO FINISH AND DO CLEAN UP
    for i=1:nConfigs
%         if submitForReal
            all_tasks_completed = SJ{eidx}(i).waitForTasksToFinish();
            if all_tasks_completed
                SJ{eidx}(i).copyBackResults();
                SJ{eidx}(i).postprocessBOTMSorting();
            end
%         else
%              SJ{eidx}(i).postprocessBOTMSorting();
%         end
    end
end

toc(t_start_all)
disp('Done Grid-Sorting!');






%% CODE SNIPETS FOR DEBUGGING
            if 0
                %% TEST PREPROCESSING
                i = 1;
                k = 1;
                cmea = mysort.mea.CMOSMEA(pathdefs.h5FileLists{i}{k});
                cpx = double(hdf5read(pathdefs.h5FileLists{i}{k}, ['/Sessions/Session0/channel_posx']));

                mysort.plot.figure([700 700]);
                n = 17.5;
                x = cmea.MultiElectrode.electrodePositions(:,1)/n;
                y = cmea.MultiElectrode.electrodePositions(:,2)/n;
                plot(x,y, 'x')

                %% IF YOU WANT TO RUN ONE TASK ON THE LOCAL MACHINE
                eidx = 1;
                i = 2;
                task = 12;
                SJ{eidx}(i).runTaskLocally(task)
            
                %% START ONLY ONE CONFIGURATION THAT RAN ALREADY ONCE (AND THE SJ{eidx} STILL EXISTS!)
                eidx = 1;
                i = 1;
                include_copy_data = true;
                if include_copy_data
                    SJ{eidx}(i).copyDataFilesToBinary(); % there is an alternative here, copyDataFiles (without "ToBinary")
                end
                
                SJ{eidx}(i).prepareTasks();
                project_token = fullfile(token_path, ['start_' PATHDEFS(eidx).jobname '.tok']);
                mysort.util.logToFile(project_token, 'Submit!');
                disp('Submitting Job to Grid-engine...');
                pause(10)
                assert(~exist(project_token, 'file'), 'Token is still there, deamon dead?!');
                disp('Submit sucessful.');

                %% WAIT ONLY FOR ONE CONFIGURATION
                eidx = 1;
                i = 1;
                expName = expNames{eidx};    
                PATHDEFS(eidx).preprocessedOutPath = fullfile(pathdefs.inputPath, [expName outSuffix]);                
                all_tasks_completed = SJ{eidx}(i).waitForTasksToFinish();
                if all_tasks_completed
                    SJ{eidx}(i).copyBackResults(pathdefs.sortOutPaths{i});
                    SJ{eidx}(i).postprocessBOTMSorting();
                end
                %% ONLY RE-CREATE job.sh files
                for i=1:length(pathdefs.h5FileLists)
                    SJ{eidx}(i).constructJobSH();
                end
                

            end