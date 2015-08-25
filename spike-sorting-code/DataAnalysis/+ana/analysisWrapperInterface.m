classdef analysisWrapperInterface < handle
    % The idea of this interface is to take away the administrative tasks
    % like folder creation, error checking, adding of trials etc from
    % analysis that are carried out on a data set which might get more and
    % more data files over time. Files are easily added by just running the
    % analysis on the files name. The derived classes must implement the
    % runTrial_(self, name) function and here, the name is the only
    % identifier that should be used to load and save data.
    properties
%         parentAnas
        analysisName
        dataInputPath
        dataOutPutPath
        dataOutPrefix
        figurePath
        figureOutPrefix
        oldFigureHandles
        
        trialInformationFile
        
        trials
        
        resultBuffer
        resultBufferName
        
        dispFun
    end
    methods(Abstract, Access=public)
        validNames = getValidTrialNames(self);
        [figureHandles figureNames] = makeFigures_(self, name, result, P);
    end
    methods(Abstract, Access=protected)
        result = runTrial_(self, name);
    end
    methods
        %------------------------------------------------------------------
        function self = analysisWrapperInterface(analysisName, dataInputPath, dataOutPutPath, figurePath)
            self.analysisName = analysisName;
            self.dataInputPath = dataInputPath;
            self.dataOutPutPath = dataOutPutPath;
            self.figurePath = figurePath;
            self.dispFun = @(x) disp([analysisName ': ' x]);
            if ~exist(self.dataOutPutPath, 'file')
                mkdir(self.dataOutPutPath);
            end
            if ~exist(self.figurePath, 'file')
                mkdir(self.figurePath);
            end
            self.dataOutPrefix = fullfile(self.dataOutPutPath, [self.analysisName '_']);
            self.figureOutPrefix = fullfile(self.figurePath, [self.analysisName '_']);
            self.trialInformationFile = [self.dataOutPrefix 'anaTIF.mat'];
            self.loadTrialInformation();
        end
        %------------------------------------------------------------------
        function result = run(self, name)
            self.loadTrialInformation();
            idx = self.getTrialIndex4TrialName(name);
            if idx == 0
                assert(self.trialNameIsValid(name), 'This is not a valid trial name! (%s)', name);
                idx = self.addTrial(name);
            end
            result = self.runIdx(idx, name);
        end 
        %------------------------------------------------------------------
        function runAllValidTrials(self)
            VN = self.getValidTrialNames();
            for i=1:length(VN)
                self.run(VN{i});
                try
                    close(self.oldFigureHandles);
                catch
                end
            end
        end
        %------------------------------------------------------------------
        function bIsValid = trialNameIsValid(self, name)
            validNames = self.getValidTrialNames();
            bIsValid = any(cellfun(@(x) strcmp(x, name), validNames));
        end
        %------------------------------------------------------------------
        function loadTrial(self, name)
            self.loadTrialInformation();
            idx = self.getTrialIndex4TrialName(name);
            if idx == 0
                self.run(name);
                idx = self.getTrialIndex4TrialName(name);
                assert(idx>0, 'Something went wrong!')
            end
            if ~(strcmp(self.trials.processState{idx}, 'done') || ...
                 strcmp(self.trials.processState{idx}, 'creatingFigures'))
                self.dispFun('Cannot load trial since it was not properly processed yet. Trying to process...')
                self.runIdx(idx, name);
            end 
        end
        %------------------------------------------------------------------
        function makeAllFigures(self, varargin)
            VN = self.getValidTrialNames();
            for i=1:length(VN)
                self.makeFigures(VN{i}, varargin{:});
            end            
        end
        %------------------------------------------------------------------
        function fh = makeFigures(self, name, varargin)
            P.saveFigArgs = {'fig', 0, 'png', 1};
            P.saveFigurePath = [];
            P.makeFigureWindows = 1;
            P = mysort.util.parseInputs(P, varargin);
            try
                close(self.oldFigureHandles);
            catch
            end
            idx = self.getTrialIndex4TrialName(name);
            if idx==0
                self.loadTrial(name);
                idx = self.getTrialIndex4TrialName(name);
                assert(idx>0, 'This is not a valid trial!');
            end
            if ~(strcmp(self.trials.processState{idx}, 'done') || ...
                 strcmp(self.trials.processState{idx}, 'creatingFigures'))
                 self.dispFun(['Trial not processed yet! (state: ' self.trials.processState{idx} ', msg: ' self.trials.failInformation{idx,1} ')']);
                 return
            end
            [figureHandles figureNames] = self.makeFigures_(name, P);
            for i=1:length(figureHandles)
                mysort.plot.savefig(figureHandles(i), ...
                    [self.figureOutPrefix '_' name '_' figureNames{i}], P.saveFigArgs{:}); 
            end
            self.oldFigureHandles = figureHandles;
        end
        %------------------------------------------------------------------
        function resetTrial(self, name)
            self.loadTrialInformation();
            idx = self.getTrialIndex4TrialName(name);
            if idx == 0
                assert(self.trialNameIsValid(name), 'Not a valid Trial name!');
                return
            end
            self.setTrialNew(idx);
            if length(self.trials.persistentVariables) < idx
                return
            end
            for i=1:length(self.trials.persistentVariables{idx})
                bfile = self.bufferFileForTrialVariable(name, self.trials.persistentVariables{idx}{i});
                if exist(bfile, 'file')
                    delete(bfile);
                end
            end
        end
        %------------------------------------------------------------------
        function resetAllTrials(self)
            VN = self.getValidTrialNames();
            for i=1:length(VN)
                self.resetTrial(VN{i});
            end           
        end    
        %------------------------------------------------------------------
        function resetAllInProcess(self)
            VN = self.getValidTrialNames();
            for i=1:length(VN)
                idx = self.getTrialIndex4TrialName(VN{i});
                if idx > 0 && strcmp(self.trials.processState{idx}, 'inProcess')
                    self.resetTrial(VN{i});
                end
            end           
        end           
        
        %------------------------------------------------------------------
        function v = getVariable(self, name, varName)
            self.loadTrial(name);
            if self.variableIsInMemory(name, varName)
                v = self.resultBuffer.(varName);
            else
                v = self.loadVariable(name, varName);
            end
        end
        %------------------------------------------------------------------
        function R = getAllVariables(self, name)
            R = [];
            self.loadTrial(name);
            idx = self.getTrialIndex4TrialName(name);
            for i=1:length(self.trials.persistentVariables{idx})
                vn = self.trials.persistentVariables{idx}{i};
                R.(vn) = self.getVariable(name, vn);
            end
        end        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    methods(Access=protected)
        %------------------------------------------------------------------
        function bfile = bufferFileForTrialVariable(self, name, varName)
            bfile = fullfile(self.dataOutPutPath, [self.analysisName '_' name '_' varName '_buffer.mat']);
        end        
        %------------------------------------------------------------------
        function b = variableIsInMemory(self, name, varName)
            b=0;
            if isfield(self.resultBufferName, varName) && strcmp(self.resultBufferName.(varName), name)
                b=1;
            end
        end
        %------------------------------------------------------------------
        function v = loadVariable(self, name, varName)
            bfile = self.bufferFileForTrialVariable(name, varName);
            assert(exist(bfile, 'file')>0, 'This variable was not computed during run!');
            d = load(bfile);
            v = d.(varName);
            self.resultBufferName.(varName) = name;
            self.resultBuffer.(varName) = v;
        end
        
%         %------------------------------------------------------------------
%         function registerParentAnalysis(self, ana)
%             assert(self.parentNotRegisteredYet(ana), 'This analysis was already registered!')
%             self.parentAnas(end+1) = ana;
%         end
%         %------------------------------------------------------------------
%         function b = parentNotRegisteredYet(self, ana)
%             b = self.parentAnaIdx4parentAnaName(ana.analysisName)>0;
%         end
%         %------------------------------------------------------------------
%         function idx = parentAnaIdx4parentAnaName(self, anaName)
%             idx = 0;
%             for i=1:length(self.parentAnas)
%                 if strcmp(self.parentAnas(i).analysisName, anaName)
%                     idx=i;
%                     return
%                 end
%             end
%         end        
        %------------------------------------------------------------------
        function idx = addTrial(self, name)
            idx = length(self.trials.names)+1;
            self.trials.names{idx} = name;
            self.trials.processState{idx} = 'new';
            self.trials.failInformation{idx} = [];
            self.trials.persistentVariables{idx} = [];
            self.saveTrialInformation(); 
        end
        %------------------------------------------------------------------
        function result = runIdx(self, idx, name)
            result = [];
            self.dispFun(sprintf('Running %s ...', name))
            if strcmp(self.trials.processState{idx}, 'done')
                self.dispFun(' already done.')
                return
            elseif strcmp(self.trials.processState{idx}, 'inProcess')
                self.dispFun(' is already in process.')
                return
            end
            
            self.setTrialInProcess(idx);
%             try
                result = self.runTrial_(name);
                if ischar(result)
                    self.runFailed(name, result);
                    result = [];
                    return
                elseif isempty(result)
                    self.runFailed(name, 'Result was empty');
                    return                    
                end
                fnames = fieldnames(result);
                self.trials.persistentVariables{idx} = fnames;
                for i=1:length(fnames)
                    self.resultBufferName.(fnames{i}) = name;
                    bfile = self.bufferFileForTrialVariable(name, fnames{i});
                    eval(sprintf('%s = result.%s;', fnames{i}, fnames{i}));
                    save(bfile, fnames{i});
                end
                self.resultBuffer = result;
                self.setTrialCreatingFigures(idx);
                self.makeFigures(name, result);
                self.setTrialDone(idx);
%             catch
%                 result = [];
%                 disp('Trial Processing failed!')
%                 errStr = mysort.util.buildLastErrString();
%                 disp(errStr);
%                 self.setTrialFailed(idx, errStr);
%             end
        end
        %------------------------------------------------------------------
        function runFailed(self, name, str)
            if nargin == 2
                str = 'none';
            end
            idx = self.getTrialIndex4TrialName(name);
            errStr = ['Run Failed in runTrial_! Message: ' str];
            self.setTrialFailed(idx, errStr)
            self.dispFun(errStr);
        end        
        %------------------------------------------------------------------
        function setTrialInProcess(self, idx)
            self.trials.processState{idx,1} = 'inProcess';
            self.saveTrialInformation();         
        end
        %------------------------------------------------------------------
        function setTrialCreatingFigures(self, idx)
            self.trials.processState{idx,1} = 'creatingFigures';
            self.trials.failInformation{idx,1} = [];
            self.saveTrialInformation();         
        end         
        %------------------------------------------------------------------
        function setTrialFailed(self, idx, errStr)
            self.trials.processState{idx,1} = 'failed';
            self.trials.failInformation{idx,1} = errStr;
            self.saveTrialInformation();         
        end 
        %------------------------------------------------------------------
        function setTrialDone(self, idx)
            self.trials.processState{idx,1} = 'done';
            self.saveTrialInformation();         
        end    
        %------------------------------------------------------------------
        function setTrialNew(self, idx)
            self.trials.processState{idx,1} = 'new';
            self.saveTrialInformation();         
        end    
        %------------------------------------------------------------------
        function loadTrialInformation(self)
            if exist(self.trialInformationFile, 'file')
                S = load(self.trialInformationFile);
                self.trials = S.S;
                if ~isfield(self.trials, 'persistentVariables')
                    self.trials.persistentVariables = cell(length(self.trials.names), 1);
                end
                return
            end         
            self.trials.names = {};
            self.trials.processState = {};
            self.trials.failInformation = {};
            self.trials.persistentVariables = {};
        end
        %------------------------------------------------------------------
        function saveTrialInformation(self)
            S = self.trials;
            save(self.trialInformationFile, 'S');         
        end 
        %------------------------------------------------------------------
        function idx = getTrialIndex4TrialName(self, names)
            if isempty(self.trials.names)
                idx= 0;
                return
            end
            if ischar(names)
                names = {names};
            end
            idx = zeros(1, length(names));
            for i=1:length(names)
                match = find(cellfun(@(x) strcmp(x, names{i}), self.trials.names),1);
                if isempty(match)
                    idx(i) = 0;
                else
                    idx(i) = match;
                end
            end
        end  
        %------------------------------------------------------------------
        function names = getTrialName4TrialIdx(self, idx)
            names = self.trials.names{idx};
        end 
    end
end
