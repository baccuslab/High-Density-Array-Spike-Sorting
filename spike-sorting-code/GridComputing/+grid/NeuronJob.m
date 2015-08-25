classdef NeuronJob < grid.GridJob
    properties (SetAccess=private)
    end
    properties
        neuronJobP
        neuronIDs
        
        neuronMap
        %DSFull
        %MES
    end
    
    methods
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = NeuronJob(jobName, templateFile, varargin)
            self = self@grid.GridJob(jobName, varargin{:});
            self.taskType = 'NeuronJob';
            
            p = struct;
            %p.useFilter = 0;
            p.nNeurons = [];
            
            p = mysort.util.parseInputs(p, self.P_untreated, 'error');
            self.neuronJobP = p;
            
            % When class is reloaded, use this:
            %if exist(self.jobName, 'file') == 2
            %disp('Re-Create DataStructure object...');
            %self.DSFull = mysort.mea.compoundMea(self.files.data, 'useFilter', self.sortJobP.useFilter, 'name', self.jobName);%'PREFILT')
            %end
            
            self.startlocation.files.data{1} = templateFile;

            
            self.completed.loadTemplateFile = false;
            self.completed.setParameters = false;
            self.completed.copyBackResults = false;
            
            self.completed.bar = false;
        end
        
        
        % -----------------------------------------------------------------
        % Load the template files and create tasks for each detected neuron:
        function loadTemplateFile(self)%, templateFile)
            if self.completed.loadTemplateFile return; end
            
            % Save check data files:
            %self.startlocation.files.data{1} = templateFile;
            assert( exist(self.startlocation.files.data{1}, 'file') == 2, 'template file missing!' );
            [pathstr,name,ext] = fileparts(self.startlocation.files.data{1});
            self.startlocation.files.data{2} = fullfile(pathstr, 'groupFile.mat');
            assert( exist(self.startlocation.files.data{2}, 'file') == 2, 'groupFile.mat missing!' );
            self.startlocation.files.data{3} = fullfile(pathstr, 'G_struct.mat');
            assert( exist(self.startlocation.files.data{3}, 'file') == 2, 'G_struct.mat missing!' );
            
            self.scratchlocation.folders.output = fullfile(self.scratchlocation.folders.main, 'output');
            [dir_exists,mess,messid] = mkdir( self.scratchlocation.folders.output );
            assert(dir_exists, 'Output directory could not be created!');
            
            % Create task for each neuron:
            load(self.startlocation.files.data{1}) % -> T_merged, gdf_merged, localSorting, localSortingID
            
            Tf = size(T_merged,1);
            nNeurons = size(T_merged,3);
            
            maxCh = 10; % maximal number of channels
            self.neuronMap = mysort.wf.tOptimalChannelsPerUnit(T_merged, maxCh);
            self.neuronIDs = unique(gdf_merged(:,1));
            
            self.startIndex = 1;
            self.endIndex = length(self.neuronIDs);
            self.taskIDs = self.startIndex:self.endIndex;
            
            % limit number of groups if necessary:
            if ~isempty(self.neuronJobP.nNeurons)
                self.neuronJobP.nNeurons = min(self.neuronJobP.nNeurons, length(self.taskIDs) );
                self.taskIDs = self.taskIDs( 1:self.neuronJobP.nNeurons) ;
                self.endIndex = self.neuronJobP.nNeurons;
            end
            self.completed.loadTemplateFile = true;
        end
        
        % -----------------------------------------------------------------
        % Varargin: strings, that produce a struct field with the value 'true'
        function setParameters(self, varargin)
            if self.completed.setParameters return; end
            
            taskParameters = struct;
            
            for i=1:length(varargin)
                taskParameters.do.(varargin{i}) = true;
            end
            
            taskParameters.runName = self.jobName;
            taskParameters.outputPath = self.scratchlocation.folders.output;
            taskParameters.Fs = 20000;
            
            %% Create cell variable allTaskParamters:
            for ii = 1:length(self.taskIDs)
                
                self.scratchlocation.folders.neuron{ii} = fullfile( taskParameters.outputPath, ['neuron' num2str(self.neuronIDs(ii)) ]);
                [dir_exists,mess,messid] = mkdir( self.scratchlocation.folders.neuron{ii} );
                assert(dir_exists, 'Output directory could not be created!');
                
                taskParameters.neuronFolder = self.scratchlocation.folders.neuron{ii};
                taskParameters.taskID = self.taskIDs(ii);
                taskParameters.neuronID = self.neuronIDs(ii);
                self.allTaskParameters{ii} = taskParameters;
            end
            
            % FINISH
            self.completed.setParameters = true;
        end
        
        % -----------------------------------------------------------------
        % Copy results from the grid folder back to the destination folder
        function copyBackResults(self, destinationFolder, copyData)
            if self.completed.copyBackResults return; end
            
            %if nargin < 3 copyData = false; end
            
            disp('Copy back result file(s)...');
            self.destinationlocation.folders.main = destinationFolder;
            [dir_exists,mess,messid] = mkdir(self.destinationlocation.folders.main);
            assert(dir_exists, 'Destination directory could not be created!');
            
            %self.files.groups = self.copyFile(self.files.groups, destinationFolder);
            %self.folders.report = self.copyFolder(self.folders.report, destinationFolder);
            %self.files.summary = self.copyFile(self.files.summary, destinationFolder);
            %self.folders.output = self.copyFolderContent(self.folders.output, destinationFolder);
            for ii = 1:length(self.taskIDs)
                self.destinationlocation.folders.neuron{ii} = self.copyFolder(self.scratchlocation.folders.neuron{ii}, destinationFolder);
            end
            
            %if copyData
            %    disp('Copy back data files...')
            %    self.folders.data = self.copyFolder(self.folders.data, destinationFolder);
            %end
            disp('Sorting results have been copied to:')
            disp(self.destinationlocation.folders.main)
            
            self.completed.copyBackResults = true;
        end
        
    end
    
    methods(Static)
        %------------------------------------------------------------------
        %------------------------- RUN FUNCTION ---------------------------
        %------------------------------------------------------------------
        function run(taskFile, debugFlag)
            if nargin < 2
                debugFlag = false;
            end
            
            %% List of default parameters:
            taskP = struct;
            taskP.do.createGraphs = false;
            taskP.do.cutSpikes = false;
            
            %% Load taskFile:
            load(taskFile);
            
            %% Write "taskParameters" to struct "taskP" and "sortP":
            %sortP = mysort.util.mergeStructs(sortP, taskParameters.sortingParameters);
            %taskParameters = rmfield(taskParameters,'sortingParameters');
            
            taskP = mysort.util.mergeStructs(taskP, taskParameters);
            clear taskParameters;
                        
            % (Re-)Set the reporting file:
            rep = mysort.ds.binaryFileMatrix(taskP.reportFile, [1 2], 'writable', true);
            rep(:,:) = [0 0];
            
            % Check that data files exist:
            for f = 1:length(taskP.dataFiles)
                assert( exist(taskP.dataFiles{f}, 'file') == 2, ['Task aborted: dataFile ' taskP.dataFiles{f} ' not found!']);
                load( taskP.dataFiles{f} )
            end
            
            try
                disp(['Starting neuron with neuronID: ' num2str( taskP.neuronID) ' ...' ])
                
                completed = struct;
                Fs = 20000
                %neuronFolder = fullfile( taskP.outputPath, ['neuron' num2str(taskP.neuronID) ]);
                %[dir_exists,mess,messid] = mkdir( neuronFolder );
                %assert(dir_exists, 'Output directory could not be created!');
                
                %--- Functions to implement: ---
                % Create graphs for each neuron
                if (taskP.do.createGraphs);
                    disp('Start creating graphs...')
                    
                    figures = struct;
                    figures.f1 = fullfile(taskP.neuronFolder, 'figure1.fig');
                    if ~exist(figures.f1, 'file')
                        disp('creating figure1')
                        F = figure;
                        ah = subplot(2,3,1)
                        p = ana.roland.rasterplot( gdf_merged(gdf_merged(:,1) == taskP.neuronID,2), 'ah', ah )
                        subplot(2,3,2)
                        plot(1:10, 1:10);
                        %mysort.plot.spiketrain(gdf_merged, 'restrictTo', taskP.neuronID)
                        subplot(2,3,3)
                        fr = ana.roland.firing_rate( gdf_merged(gdf_merged(1:1000,1) == taskP.neuronID,2)/taskP.Fs, 'sigma', 0.1)
                        plot(fr.time, fr.firing_rate{1})
                        xlabel('time [s]'); ylabel('spiking rate')
                        
                        subplot(2,3,4)
                        mysort.plot.isiViolations(0, gdf_merged, 3.0*Fs)
                        
                        mysort.plot.savefig(F, figures.f1)
                    end
                    % -> function cutSpikes(...)
                end
                
                
                % Cut spikes on the grid for each neuron and produce an H5 file
                if (taskP.do.cutSpikes)
                    disp('Start cutting spikes...')
                    % -> function createGraphs(...)
                end
                
                if (taskP.do.test);
                    disp('Start test ...')
                    
                    x0 = taskP.neuronID;
                    save( fullfile(taskP.neuronFolder, 'test.mat'), 'x0')
                end
                
                
                %% Write to reporter file:
                disp('Writing results...')
                rep = mysort.ds.binaryFileMatrix(taskP.reportFile, [1 2], 'writable', true);
                rep(:,:) = [1 0];
                
            catch
                disp('Catch error...')
                errStr = mysort.util.buildLastErrString();
                disp(errStr)
                
                %% Report error event to reporter file:
                rep = mysort.ds.binaryFileMatrix(taskP.reportFile, [1 2], 'writable', true);
                rep(:,:) = [0 1];
            end
        end
        
    end
end
