% Things necessary before you can run this script:
% TO DO ONCE:
% A) Become a Grid User -> talk to IT
% B) Make a checkout of the matlab SVN folder on your new Grid Home
% B-1) For this you will need to set up your SSH key on the submit host to
%      be able to access the SVN
% C) copy the startup.m from the matlab folder to your Grid Home and set
%    the correct paths
% D) Test this by logging into the submit system, loading the matlab
%    module, starting matlab and see if the mysort package is loaded
% E) ...
%
%
% TO DO EACH TIME YOU RUN A NEW JOB
% A) ...
classdef GridJob < handle
    properties (SetAccess=private)
    end
    properties
        gc
        jobName
        allTaskParameters
        %folders
        logFolder
        %files
        startlocation
        scratchlocation
        destinationlocation
        
        m_nTasksCompletedLastTime
        m_tTimeLastChange
        m_tTimeStart
        completed
        taskIDs
        taskType
        
        startIndex
        endIndex
        
        P
        P_untreated
    end
    
    methods
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = GridJob(job_name, varargin)
            % Starts a grid job by running a number of parallel matlab
            % calls on a series of config files. Each individual matlab task
            % will load a config file which is unique for this task and run
            % the same script as all other tasks as well.
            self.gc = grid.config();
            
            is_letter = isstrprop(job_name, 'alpha');
            assert(is_letter(1), 'Job name must begin with a letter!')
            
            p.runLocallyOn = [];
            p.queue = 'regular';
            p.resume = false;
            p.dataPath = [];
            p.logFolder = ['~/log/' job_name '/'];
            
            [p p2] = mysort.util.parseInputs(p, varargin, 'split');
            self.P_untreated = p2; self.P = p;
            
            % For debugging on a local machine:
            if ~isempty(self.P.runLocallyOn) self.gc.home = self.P.runLocallyOn; end
            
            
            self.jobName = job_name;
            self.logFolder = p.logFolder;
            self.createFolders();
            
            self.completed.prepareTasks = false;
            self.completed.copyDataFilesToScratch = false;
            self.completed.createTaskFiles = false;
            self.completed.waitForTasksToFinish = false;
        end
        
        function createFolders(self)
            
            if ~isempty(self.P.dataPath)
                %% Put sorting output at the same location where the preprocessed files are
                [dir_exists,mess,messid] = mkdir(self.P.dataPath, self.jobName);
                assert(dir_exists, 'Job directory could not be created!');
                self.scratchlocation.folders.main = self.P.dataPath, self.jobName;
            else
                %% Create job folder:
                [dir_exists,mess,messid] = mkdir(self.gc.home, self.jobName);
                assert(dir_exists, 'Job directory could not be created!');
                self.scratchlocation.folders.main = fullfile(self.gc.home, self.jobName);
                
                %% Data folder:
                self.scratchlocation.folders.data = fullfile(self.scratchlocation.folders.main, 'data');
                [dir_exists,mess,messid] = mkdir(self.scratchlocation.folders.data);
                assert(dir_exists, 'Data file directory could not be created!');
            end
            %% LOG folder:
            if ~exist(self.logFolder, 'file')
                mkdir(self.logFolder);
            end
            
            %% Task files folder:
            self.scratchlocation.folders.tasks = fullfile( self.scratchlocation.folders.main, 'taskFiles');
            [dir_exists,mess,messid] = mkdir( self.scratchlocation.folders.tasks);
            assert(dir_exists, 'Task file directory could not be created!');
            
            %% Report files folder:
            self.scratchlocation.folders.report = fullfile( self.scratchlocation.folders.main, 'reportFiles');
            [dir_exists,mess,messid] = mkdir(self.scratchlocation.folders.report);
            assert(dir_exists, 'Report file directory could not be created!');
        end
        
        % -----------------------------------------------------------------
        function prepareTasks(self)%, varargin)
            if self.completed.prepareTasks return; end
            %self.P = mysort.util.parseInputs(self.P, varargin, 'error');
            
            %self.copyDataFilesToScratch();
            self.createReportFiles();
            self.createTaskFiles();
            self.constructJobSH();
            
            % Set resume value to true
            self.setResume(true);
            self.completed.prepareTasks = true;
        end
        
        % -----------------------------------------------------------------
        function setResume(self, r)
            self.P.resume = r;
        end
        
        % -----------------------------------------------------------------
        function fileName = saveToFile(self)
            fileName = fullfile(self.folders.main,'GridJob.mat');
            save(fileName,'self');
        end
        
        % -----------------------------------------------------------------
        % This function creates a file for each task, containing only one
        % structure named taskParameters. The class element
        % allTaskParameters{:} has to be specified before.
        function createTaskFiles(self)
            if self.completed.createTaskFiles return; end
            
            self.scratchlocation.files.tasks = fullfile( self.scratchlocation.folders.tasks, '/taskFile');
            
            assert( self.nTasks == length(self.allTaskParameters), 'Error: number of tasks does not correspond to size of allTaskParameters!');
            for ii = 1:self.nTasks
                taskParameters = self.allTaskParameters{ii};
                taskParameters.reportFile = self.scratchlocation.files.report{ii};
                taskParameters.reportFolder = self.scratchlocation.folders.report;
                taskParameters.dataFiles = self.scratchlocation.files.data;
                
                taskType = self.taskType;
                save( [self.scratchlocation.files.tasks num2str(taskParameters.taskID)], 'taskParameters', 'taskType');
            end
            self.completed.createTaskFiles = true;
        end
        
        % -----------------------------------------------------------------
        function createReportFiles(self)
            for ii = 1:self.nTasks
                self.scratchlocation.files.report{ii} = fullfile( self.scratchlocation.folders.report, ['task' num2str(self.taskIDs(ii)) '.dat']);
                %
                %                 if  ~self.P.resume
                %                     if exist(self.scratchlocation.files.report{ii}, 'file') == 2
                %                         delete( self.scratchlocation.files.report{ii} );
                %                     end
                %
                %                     rep = mysort.ds.binaryFileMatrix(self.scratchlocation.files.report{ii}, [1 2], 'writable', true);
                %                     rep(:,:) = [0 0];
                %                 end
                
                % Only create new report files when it doesn't exist yet:
                if exist(self.scratchlocation.files.report{ii}, 'file') ~= 2
                    rep = mysort.ds.binaryFileMatrix(self.scratchlocation.files.report{ii}, [1 2], 'writable', true);
                    rep(:,:) = [0 0];
                end
            end
        end
        
        % -----------------------------------------------------------------
        function out = waitForTasksToFinish(self, pause_duration)
            %Wait for all tasks to finish
            % OUT = waitForTasksToFinish(self, PAUSE_DURATION=10)
            % This function regularly checks if the tasks have finished and
            % displays the number of successfully finished tasks and the
            % number of errors.
            %  - returns 'true' when all tasks have completed.
            %  - returns 'false' when there has been at least one error.
            %
            %  PAUSE_DURATION: waiting time between loops in seconds.
            
            if self.completed.waitForTasksToFinish out = true; return; end
            
            if nargin < 2 pause_duration = 10; end
            disp([self.jobName ': Waiting for tasks to finish (wait for ' num2str(pause_duration) ' seconds)...']);
            
            while true
                nCompleted = 0;
                nErrors = 0;
                for ii = 1:self.nTasks
                    
                    try
                        rep = mysort.ds.binaryFileMatrix(self.scratchlocation.files.report{ii}, [1 2], 'writable', false);
                    catch
                        disp([self.jobName ': Binary file for task ' num2str( self.taskIDs(ii) ) ' threw exception!']);
                        rep = [0 1];
                    end
                    
                    if rep(1,2) > 0
                        disp([self.jobName ': Error in task ' num2str( self.taskIDs(ii) ) ]);
                        nErrors = nErrors + 1;
                    end
                    if rep(1,1) > 0
                        nCompleted = nCompleted + 1;
                    end
                end
                diff_completed = nCompleted - self.m_nTasksCompletedLastTime;
                self.m_nTasksCompletedLastTime = nCompleted;
                try
                    dt = toc(self.m_tTimeLastChange);
                catch
                    dt = 0;
                end
                minutes = floor(dt/60);
                sec = round(dt-minutes*60);
                dtstr = sprintf('%d min %d sec', minutes, sec);
                if diff_completed>0
                    self.m_tTimeLastChange = tic;
                end
                fprintf('%s: %d (+ %d) out of %d tasks completed. (Time to last Change: %s)\n', self.jobName, nCompleted, diff_completed, self.nTasks, dtstr);
                
                if nCompleted + nErrors == self.nTasks break; end
                pause(pause_duration);
            end
            
            out = true;
            if ~nErrors
                disp([self.jobName ': All tasks completed successfully...']);
            else
                disp([self.jobName ': ATTENTION: Waiting ended without completing all tasks!']);
                disp([num2str(nErrors) ' errors detected!']);
                out = false;
            end
            self.summarizeReports();
            
            self.completed.waitForTasksToFinish = true;
        end
        
        % -----------------------------------------------------------------
        function [bFinished, nErrors, nUnfinished, nCompleted, nTotal] = isTaskFinished(self)
            % Checks the status of all tasks 
            nTotal = self.nTasks;
            bFinished = true;
            nCompleted = 0;
            nErrors = 0;
            
            for ii = 1:self.nTasks
                try
                    rep = mysort.ds.binaryFileMatrix(self.scratchlocation.files.report{ii}, [1 2], 'writable', false);
                catch
                    disp([self.jobName ': Binary file for task ' num2str( self.taskIDs(ii) ) ' threw exception!']);
                    rep = [0 1];
                end

                if rep(1,2) > 0
                    disp([self.jobName ': Error in task ' num2str( self.taskIDs(ii) ) ]);
                    nErrors = nErrors + 1;
                end
                if rep(1,1) > 0
                    nCompleted = nCompleted + 1;
                end
            end
            
            nUnfinished = nTotal - nErrors - nCompleted;
            
            if nUnfinished > 0
                bFinished = false;
            end
        end        
        
        % -----------------------------------------------------------------
        function summarizeReports(self)
            % This function reads all reports files and creates a summary of
            % them.
            
            %% Check for most recent job id:
            try
                currentJobId = fullfile(self.scratchlocation.folders.report, 'currentJobID.txt')
                fid = fopen(currentJobId);
                job_id_str = fread(fid, '*char')'; fclose(fid);
            catch
                job_id_str = 'undetermined_job';
            end
            
            %% Create a summary file for this job:
            self.scratchlocation.files.summary = fullfile(self.scratchlocation.folders.main, [job_id_str '_summary.txt'] );
            self.scratchlocation.files.summaryMAT = fullfile(self.scratchlocation.folders.main, [job_id_str '_summary.mat'] );
            summary = struct;
            
            fid = fopen(self.scratchlocation.files.summary, 'a+');
            fprintf(fid, ['Job: ' job_id_str '\n'])
            
            summary.tasksError = [];
            summary.nCompleted = 0;
            summary.times = {};
            for ii = 1:self.nTasks
                rep = mysort.ds.binaryFileMatrix(self.scratchlocation.files.report{ii}, [1 2], 'writable', false);
                if rep(1,2) > 0
                    summary.tasksError = [summary.tasksError self.taskIDs(ii)];
                end
                if rep(1,1) > 0
                    summary.nCompleted = summary.nCompleted + 1;
                end
                
                %% Extract the time of this task:
                try
                    startTimeFile = fullfile(self.scratchlocation.folders.report, ['startTime'  job_id_str '_'  num2str(self.taskIDs(ii)) '.mat']);
                    load(startTimeFile);
                    
                    endTimeFile = fullfile(self.scratchlocation.folders.report, ['endTime'  job_id_str '_' num2str(self.taskIDs(ii)) '.mat']);
                    load(endTimeFile);
                    
                    summary.times{ii, 1} = startTime;
                    summary.times{ii, 2} = endTime;
                    summary.times{ii, 3} = etime(endTime, startTime);
                catch
                    summary.times{ii, 1} = clock;
                    summary.times{ii, 2} = clock;
                    summary.times{ii, 3} = 0;
                end
            end
            summary.totalDuration = self.findTotalDuration( summary.times );
            fprintf(fid, ['Total Duration (in seconds): ' num2str(summary.totalDuration) '\n']);

            fprintf(fid, ['Number of tasks completed: ' num2str(summary.nCompleted) '\n']);
            fprintf(fid, ['Number of errors:          ' num2str( length(summary.tasksError) ) '\n']);
            
            if length(summary.tasksError) > 0
                fprintf(fid, ['Errors in tasks: ' num2str( summary.tasksError ) '\n']);
            end
            fclose(fid);
            
            save(self.scratchlocation.files.summaryMAT, 'summary'); 
        end
        
        function totalDuration = findTotalDuration(self, times )
            minTime = times{1,1};
            maxTime = times{1,2};
            for ii = 1:length(times)
                startT = times{ii,1};
                stopT = times{ii,2};
                
                if etime(startT, minTime) < 0
                    minTime = startT;
                end
                if etime(stopT, maxTime) > 0
                    maxTime = stopT;
                end
            end
            
            totalDuration = etime(maxTime, minTime);
        end
        
        % -----------------------------------------------------------------
        function str = constructConsoleJobStartCommand(self)
            str = sprintf('qsub -q %s.q %s_job.sh', self.P.queue, self.jobName);
        end
        
        % -----------------------------------------------------------------
        function filename = constructJobSH(self)
            
            str = '#!/bin/bash\n';
            str = [str '#$ -V\n'];
            str = [str '#$ -cwd\n'];
            str = [str '#$ -o ' self.logFolder '\n'];
            str = [str '#$ -e ' self.logFolder '\n'];
            str = [str '#$ -q ' self.P.queue '.q\n'];
            str = [str sprintf('#$ -t %d-%d\n', self.startIndex, self.endIndex)];
            str = [str sprintf('matlab -nodisplay -r "grid.GridJob.runTask(''%s\''); exit();"', self.scratchlocation.files.tasks)];
            
            %filename = fullfile( self.gc.home, self.jobName, [self.jobName '_job.sh']);
            filename = fullfile(self.scratchlocation.folders.main, [self.jobName '_job.sh']);
            fh = fopen(filename, 'w+');
            fprintf(fh, str );
            fclose(fh);
        end
        
        % -----------------------------------------------------------------
        function copyDataFilesToScratch(self)
            if self.completed.copyDataFilesToScratch return; end
            
            disp('Copy data file(s)...');
            for i = 1:length(self.startlocation.files.data)
                self.scratchlocation.files.data{i} = self.copyFile(self.startlocation.files.data{i}, self.scratchlocation.folders.data);
                disp([ num2str(i) ' of ' num2str(length(self.startlocation.files.data)) ' copied.']);
            end
            
            self.completed.copyDataFilesToScratch = true;
        end
        
        
        % -----------------------------------------------------------------
        % Copy one file to destination folder and return new file name:
        function fileName = copyFile(self, fileName, destFolder, forceCopy)
            if nargin < 4 forceCopy = false; end
            flags = '-q';%  '-vP'
            
            if forceCopy
                system(['rsync ' flags ' ' fileName ' ' destFolder]);
            else
                system(['rsync ' flags ' --ignore-existing ' fileName ' ' destFolder]);
            end
            
            [pathstr,name,ext] = fileparts(fileName);
            fileName = fullfile(destFolder, [name ext]);
        end
        
        function folderName = copyFolder(self, folderName, destFolder, forceCopy)
            if nargin < 4 forceCopy = false; end
            
            if forceCopy
                system(['rsync -vP -r' folderName ' ' destFolder]);
            else
                system(['rsync -vP -r --ignore-existing ' folderName ' ' destFolder]);
            end
            
            folderName = fullfile(destFolder, folderName);
        end
        
        function folderName = copyFolderContent(self, folderName, destFolder, forceCopy)
            if nargin < 4 forceCopy = false; end
            
            if forceCopy
                system(['rsync -vP -r' folderName '/* ' destFolder]);
            else
                system(['rsync -vP -r --ignore-existing ' folderName '/* ' destFolder]);
            end
            
            folderName = fullfile(destFolder);
        end
        % -----------------------------------------------------------------
        
        
        
        %         % Copy a folder with its contents. Optimized for speed, i.e. not
        %         % re-copying already existing files.
        %         function folderName = copyFolder(self, folderName, destFolder, contentOnly)
        %             if nargin < 4 contentOnly = false; end
        %             disp([ 'copyFolder( ' folderName ', ' destFolder ', ' num2str(contentOnly) ' )'])
        %
        %             if ~contentOnly
        %                 [pathstr,name] = fileparts(folderName);
        %                 destFolder = fullfile(destFolder, name)
        %                 if exist(destFolder, 'dir') ~= 7
        %                     [dir_exists,mess,messid] = mkdir(destFolder);
        %                 end
        %             end
        %
        %             list = dir(folderName);
        %             for i = 1:length(list)
        %                 if strcmp( list(i).name(1), '.') continue; end
        %                 f = fullfile(folderName,list(i).name);
        %
        %                 if list(i).isdir
        %                     self.copyFolder(f, destFolder);
        %                 else
        %                     self.copyFile(f, destFolder);
        %                 end
        %             end
        %         end
        %
        %         function folderName = copyFolderContent(self, folderName, destFolder)
        %             folderName = self.copyFolder(folderName, destFolder, true);
        %         end
        
        % -----------------------------------------------------------------
        function N = nTasks(self)
            N = self.endIndex - self.startIndex + 1;
        end
        
        
        % -----------------------------------------------------------------
        %              FOR DEBUGGING:
        % -----------------------------------------------------------------
        function runTaskLocally(self, task_id)
            
            setenv('SGE_TASK_ID', num2str(task_id))
            
            
            setenv('JOB_ID', '0')
            jobFolder = fullfile(self.scratchlocation.folders.report, '0')
            try
                rmdir(jobFolder, 's');
            catch
            end
            mkdir(jobFolder);
            grid.GridJob.runTask(self.scratchlocation.files.tasks, true);
        end
        
        
        function prepareTest(self)
            self.startIndex = 3;
            self.endIndex = 12;
            %self.nTasks = self.endIndex - self.startIndex + 1;
            
            self.files.tasks = fullfile( self.gc.home, self.jobName, '/taskFile');
            
            var1 = 'var1';
            var2 = 'var2';
            num1 = 12345;
            nTasks = self.numTasks();
            outputPath = fullfile( self.gc.home, self.jobName, 'results');
            
            %             %% Create folder, then task files
            %             [dir_exists,mess,messid] = mkdir(self.gc.home, self.jobName);
            %             assert(dir_exists, 'Job Directory could not be created!');
            
            for task_id = self.startIndex:self.endIndex
                save( [self.files.tasks num2str(task_id)], 'outputPath', 'var1', 'var2', 'num1', 'task_id', 'nTasks');
            end
            
            test = 1;
            self.constructJobSH(test);
        end
    end
    
    methods (Abstract, Static)
        run(taskFile)
    end
    
    methods (Static)
        
        function out = runTask(taskName, debugFlag)
            %Static function of GridJob called by each node on the grid.
            % OUT = grid.GridJob.runTask(TASKNAME, DEBUGFLAG)
            %
            % OUT is boolean indicating success or failure of task.
            disp(['Starting grid.GridJob.runTask(' taskName ') ...'])
            
            if nargin < 2
                debugFlag = false;
            end
            
            %% Find task-ID and load its taskFile:
            task_id_str = getenv('SGE_TASK_ID');
            assert( ~isempty(task_id_str), 'No SGE_TASK_ID found!');
            taskFile = fullfile([taskName task_id_str '.mat']);
            f = load(taskFile);
            
            %% Check if the same task in the same job already exists:
            job_id_str = getenv('JOB_ID');
            processFolder = fullfile( f.taskParameters.reportFolder, job_id_str)
            [dir_exists,mess,messid] = mkdir(processFolder);
            uniqueProcessFile = fullfile( processFolder, ['task' task_id_str 'exists.dat'])
            if exist(uniqueProcessFile, 'file') == 2
                error(['Task ' task_id_str ' in job ' job_id_str ' already exists!'])
            end
            fileID = fopen(uniqueProcessFile,'w+'); fwrite(fileID, zeros(1), 'int8'); fclose(fileID);
            disp('uniqueProcessFile created.')
            
            %% Create a file that stores the current job_id:
            fileID = fopen( fullfile(f.taskParameters.reportFolder, 'currentJobID.txt'),'w+'); fwrite(fileID, job_id_str); fclose(fileID);
            clear fileID
            disp('currentJobID.txt file created.')
            
            %% Create a file that stores the start time vector:
            startTime = clock;
            save(  fullfile(f.taskParameters.reportFolder, ['startTime' job_id_str '_' task_id_str '.mat']), 'startTime' );
            disp(['startTime' task_id_str '.mat file created.'])
            
            try
                grid.(f.taskType).run(taskFile, debugFlag)
                out = true;
            catch
                disp( ['Possible error source: No taskType found: ' f.taskType ]);
                out = false;
                return;
            end
            
            %% Create a file that stores the current time vector:
            endTime = clock;
            save(  fullfile(f.taskParameters.reportFolder, ['endTime'  job_id_str '_'  task_id_str '.mat']), 'endTime' );
            disp(['endTime' task_id_str '.mat file created.'])
            
            
            disp('grid.GridJob.runTask() finished.')
        end
        
    end
    
end


%%
% % Initialize Console:
% module load repo/grid
% module load grid/grid
% module load matlab/R2014a

% % .sh script to actually run all the jobs
% #!/bin/bash
% #$ -cwd
% #$ -t 1-50
% time echo "scale=5000; 4*a(1)" | bc -l

% % console command to trigger the job. "regular.q" is defining the queue
% % in which to run the job
% qsub -q regular.q job.sh
