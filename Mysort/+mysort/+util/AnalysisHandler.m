
classdef AnalysisHandler < mysort.util.DebuggableClass
    properties (Constant)

    end
    properties
        name
        bLoaded
        eval_function
        save_path
        parent
        
        fullpath
        bufferFile
        resultsInSingleFiles
        
        parameterCounts
        combinatorics

        saveItems
    end
    
    methods (Abstract)

    end
    
    methods
        %%% ------------------------------------------------------
        function self = AnalysisHandler(name, eval_function, parameterNames,...
                parameterLists, savePathOrParent, varargin)
            self = self@mysort.util.DebuggableClass(varargin{:});
            self.P.version = [];
            self.P.resultsInSingleFiles = 0;
            self.P.keepResultsInMemory = 0;
            self.P.ignoreParameterForLoadControl = [];
            self.P = mysort.util.parseInputs(self.P, 'AnalyseHandler', varargin);
            
            self.name = name; self.P.debug_name = name;
            self.bLoaded = 0;
            self.resultsInSingleFiles = self.P.resultsInSingleFiles;
            self.eval_function  = eval_function;
            self.saveItems.parameterNames = parameterNames;
            self.saveItems.parameterLists = parameterLists;
            self.parent = [];
            if isa(savePathOrParent, 'mysort.util.AnalysisHandler')
                self.parent = savePathOrParent;
                self.save_path = self.parent.save_path;
            else
                self.save_path = savePathOrParent;            
            end
            self.fullpath = fullfile(self.save_path, self.name);
            self.bufferFile = fullfile(self.fullpath, 'buffer.mat');
            self.parameterCounts = cellfun(@length, self.saveItems.parameterLists);
            
            self.combinatorics = mysort.util.buildCombinatorics(self.getParameterCounts());
            self.saveItems.exceptionList = [];
            
            if ~exist(self.fullpath, 'dir')
                self.debugout('AnalysePath not found, creating...',self.LEVEL_EXTREME);
                mkdir(self.fullpath);
            else
                self.debugout('AnalysePath found, loading...',self.LEVEL_EXTREME);
                self.load();
            end
            
            if ~self.bLoaded
                self.saveItems.RES = {};
                self.saveItems.exceptionList = zeros(size(self.combinatorics,1),1);
            end
        end
        
        %%% ------------------------------------------------------
        function start(self)
            self.debugout('Starting...',self.LEVEL_EXTREME);
            nRES = length(self.saveItems.RES);
            nCombis       = size(self.combinatorics,1);
            nMyParameters = length(self.saveItems.parameterLists);
            timer = mysort.util.ProcessTimer(nCombis);
            for i = 1:nCombis
                timer.next();
                str = timer.getProgressString();
                self.debugout(str,self.LEVEL_EXTREME);
                calculateThis = 0;
                if self.saveItems.exceptionList(i) == 1
                    self.debugout('Exception found.', self.LEVEL_EXTREME);
                    calculateThis = 1;
                elseif self.resultsInSingleFiles
                    if ~self.singleFileIsLoadable(i);
                    	calculateThis = 1;
                    end
                else
                    if nRES < i || isempty(self.saveItems.RES{i})
                        calculateThis = 1;
                    end
                end
                
                if ~calculateThis 
                    self.debugout('Buffer found!', self.LEVEL_EXTREME);
                    continue
                end
                
                % Build value list for this parameter combination
                vals = self.buildValueList(i);
                self.debugout('Calculating...', self.LEVEL_EXTREME);
                % Execute the function with the corresponding value pairs
%                 try
                    savedir = fullfile(self.fullpath, sprintf('%04d', i));
                    if ~exist(savedir, 'dir')
                        mkdir(savedir);
                    end
                    self.saveItems.RES{i} = self.eval_function(vals{:}, savedir);
                    self.saveItems.exceptionList(i) = 0;
%                 catch ME
%                     self.debugout('Exception caught!', self.LEVEL_EXTREME);
%                     self.saveItems.RES{i} = ME;
%                     self.saveItems.exceptionList(i) = 1;
%                 end   
                if self.resultsInSingleFiles
                    self.saveSingleFile(i);
                    if ~self.P.keepResultsInMemory
                        self.saveItems.RES{i} = [];
                    end
                end
            end
            self.save();
        end
        
        %%% ------------------------------------------------------
        function vals = buildValueList(self, i)
            combi = self.combinatorics(i,:);
            vals = self.getParentResultList(i);
            if ~isempty(vals)
                tmp = self.getParameterListForCombination(combi(2:end));
                vals = {vals{:} tmp{:}};
            else
                vals = self.getParameterListForCombination(combi);
            end
        end
        
        %%% ------------------------------------------------------
        function results = getParentResultList(self, i)
            results = [];
            if ~isempty(self.parent)
                parentIdx = self.combinatorics(i,1);
                parentRESi = self.parent.getResult(parentIdx);
                tmp = self.parent.getParentResultList(parentIdx);
                
                if ~isempty(tmp)
                    results = {tmp{:} parentRESi};
                else
                    results = {parentRESi};
                end
            end            
        end
        
        %%% ------------------------------------------------------
        function vals = getParameterListForCombination(self, combi)
            vals = {};
            for k=1:length(self.saveItems.parameterLists)
                pList = self.saveItems.parameterLists{k};
                if iscell(pList)
                    vals = {vals{:} pList{combi(k)}};
                else
                    vals = {vals{:} pList(combi(k))};
                end
            end              
        end
        
        %%% ------------------------------------------------------
        function var = get(self, i, varname)
            var = self.getResult(i, varname);
            var = var.(varname);
        end
        %%% ------------------------------------------------------
        function RESi = getResult(self, i, varname)
            if nargin < 3
                varname = [];
            end
            
            if self.resultsInSingleFiles
                if self.P.keepResultsInMemory && ...
                        length(self.saveItems.RES)>=i ...
                        && ~isempty(self.saveItems.RES{i}) ...
                        && (isempty(varname) || isfield(self.saveItems.RES{i}, varname))
                    % Result is already in memory
                    if isempty(varname)
                        RESi = self.saveItems.RES{i};
                    else
                        RESi.(varname) = self.saveItems.RES{i}.(varname);
                    end
                elseif self.P.keepResultsInMemory
                    % Result should be in memory but is not
                    if isempty(varname)
                        self.saveItems.RES{i} = self.loadSingleFile(i);
                        RESi = self.saveItems.RES{i};
                    else
                        self.saveItems.RES{i}.(varname) = self.loadSingleFile(i, varname);
                        RESi.(varname) = self.saveItems.RES{i}.(varname);
                    end
                else
                    % Result should not stay in memory but has to be loaded
                    if isempty(varname)
                        RESi = self.loadSingleFile(i);
                    else
                        RESi.(varname) = self.loadSingleFile(i, varname);
                    end                         
                end
            else
                % Result is ready
                if isempty(varname)
                    RESi = self.saveItems.RES{i};
                else
                    RESi.(varname) = self.saveItems.RES{i}.(varname);
                end                
            end
        end
        
        %%% ------------------------------------------------------
        function load(self)
            % Check if buffered result file exists
            if ~exist(self.bufferFile,'file')
                self.debugout('...buffer file missing.',self.LEVEL_EXTREME);
                return
            end
            % Check if the file has the same version as this analysis
            tmp = load(self.bufferFile, 'parameterLists', 'version');
            if self.P.version ~= tmp.version
                self.debugout('...version mismatch.',self.LEVEL_EXTREME);
                return
            end
            
            % Check if the parameterLists are identical
            if ~self.paramterListsIdentical(tmp.parameterLists, self.saveItems.parameterLists)
                self.debugout('...parameter list mismatch.',self.LEVEL_EXTREME);
                return
            end
            
            % Load buffered file
            buffer = load(self.bufferFile);
            self.saveItems.RES = buffer.RES;
            self.saveItems.exceptionList = buffer.exceptionList;
            clear buffer
            self.bLoaded = 1;
            self.debugout('...loaded.',self.LEVEL_EXTREME);
        end

        %%% ------------------------------------------------------
        function b = paramterListsIdentical(self, list1, list2)
            b = 0;
            for i=1:length(list1)
                if any(self.P.ignoreParameterForLoadControl==i)
                    continue
                end
                if length(list1{i}) ~= length(list2{i})
                    return
                end
                for k=1:length(list1{i})
                    if iscell(list1{i})
                        if length(list1{i}{k}) ~= ...
                           length(list2{i}{k})
                            d = 0;
                        else
                            d = list1{i}{k} == list2{i}{k};
                        end
                    else
                        if length(list1{i}(k)) ~= ...
                           length(list2{i}(k))
                            d = 0;
                        elseif isstruct(list1{i}(k)) || ...
                               isstruct(list2{i}(k))
                            d = isstruct(list1{i}(k)) && ...
                                isstruct(list2{i}(k));
                        else                        
                            d = list1{i}(k) == list2{i}(k);
                        end
                    end
                    if any(d==0)
                        return
                    end
                end
            end
            b=1;
        end
        
        %%% ------------------------------------------------------
        function RESi = loadSingleFile(self, i, varname)
            if nargin < 3
                varname = [];
            end
            RESi = [];
            fname =  fullfile(self.fullpath, sprintf('buf%04d.h5',i));
            if self.singleFileIsLoadable(i)
                if ~isempty(varname)
                    RESi = mysort.h5.recursiveLoad(fname, ['/RESi/' varname]);
                else
                    RESi = mysort.h5.recursiveLoad(fname, '/RESi');
                end
            else
                error('File is not loadable! (%s)', fname);
            end
        end
        
        %%% ------------------------------------------------------
        function b = singleFileIsLoadable(self, i)
            b = 1;
            fname =  fullfile(self.fullpath, sprintf('buf%04d.h5',i));
            if ~exist(fname,'file')
                b = 0;
                return
            end    
%             version = hdf5read(fname, '/version_');
%             if ~isempty(self.P.version) && version ~= self.P.version
%                 b = 0;
%             end
        end
            
        %%% ------------------------------------------------------
        function save(self)
            self.saveItems.version = self.P.version;
            self.saveItems.date = datestr(now);
            saveItems = self.saveItems;
            if self.resultsInSingleFiles
                saveItems.RES = {};
            end
            save(self.bufferFile, '-struct', 'saveItems','-v7.3');
%             save(file, '-append', '-struct', '', '-v7.3');            
        end
        
        %%% ------------------------------------------------------
        function saveSingleFile(self, i)
            fname = fullfile(self.fullpath, sprintf('buf%04d.h5',i));
            version = self.P.version;
            date = datestr(now);
            RESi = self.saveItems.RES{i};
            [nCombis nParas] = size(self.combinatorics);
            hdf5write(fname, '/date', date);
            hdf5write(fname, '/nCombis', nCombis, 'WriteMode', 'append');
            hdf5write(fname, '/nParas', nParas, 'WriteMode', 'append');
            hdf5write(fname, '/version', version, 'WriteMode', 'append');
            mysort.h5.recursiveSave(fname, RESi, '/RESi')
        end
        %%% ------------------------------------------------------
        %%% -----------------HELPER FUNCTIONS---------------------
        %%% ------------------------------------------------------
        function [M ndimVar] = getResultMatrix(self, varname, wrapFunc, constraints)
            ndimVar = self.getResultVarOffset(varname, wrapFunc);
            C = self.getAllParentResolvedCombinatorics();
            Cidx = (1:size(self.combinatorics,1))';
            for i=1:size(constraints,1)
                constraintParameterName  = constraints(i,1);
                [constraintParameterNameIdx c_offset] = ...
                    self.parameterNamesToIndices(constraintParameterName, varname, wrapFunc);

                constraintParameterValueIdx = self.getParameterValueIdx(c_offset, constraintParameterNameIdx, constraints{i,2});
                myIdxSet = C(:,constraintParameterNameIdx-c_offset)==constraintParameterValueIdx;
                C = C(myIdxSet,:);
                Cidx = Cidx(myIdxSet);
            end              
            
            for i=1:size(Cidx,1)
                RESi = self.getResult(Cidx(i), varname);
                M(i,:) = [wrapFunc(RESi.(varname)) self.getParentResolvedCombinatorics(Cidx(i))];
            end
        end
        
        %%% ------------------------------------------------------
        function R  = getResultCell(self, constraints)
            if ~exist('constraints', 'var')
                constraints = [];
            end

            R = {};
            R.combis = self.getAllParentResolvedCombinatorics();
            R.parameter_names = self.getLabelList();
            R.parameter_lists = self.getParameterLists();
            R.combi_idx = (1:size(self.combinatorics,1))';
            keepidx = ones(size(self.combinatorics,1),1);
            for i=1:2:length(constraints)
                constraint_para_idx  = find(strcmp(R.parameter_names, constraints{i}));
                constraint_val_idx   = find(R.parameter_lists{constraint_para_idx}==constraints{i+1});
                keepidx = keepidx & (R.combis(:,constraint_para_idx) == constraint_val_idx);
            end
            R.combis = R.combis(keepidx,:);
            R.combi_idx = R.combi_idx(keepidx);
            for i=1:length(R.combi_idx)
                R.R{i} = self.getResult(R.combi_idx(i));
            end              
        end
        
        %%% ------------------------------------------------------
        function C = getAllParentResolvedCombinatorics(self)
            for i=1:size(self.combinatorics,1)
                C(i,:) = self.getParentResolvedCombinatorics(i);
            end
        end
        
        %%% ------------------------------------------------------
        function comb = getParentResolvedCombinatorics(self, i)
            if isempty(self.parent)
                comb = self.combinatorics(i,:);
            else
                parentIdx = self.combinatorics(i,1);
                comb = [self.parent.getParentResolvedCombinatorics(parentIdx)...
                        self.combinatorics(i,2:end)];
            end
        end
        
        %%% ------------------------------------------------------
        function offset = getResultVarOffset(self, varname, wrapfunc)
            RESi = self.getResult(1, varname);
            offset = length(wrapfunc(RESi.(varname)));
        end
        
        %%% ------------------------------------------------------
        function labels = getLabelList(self)
            labels = self.saveItems.parameterNames;
            if ~isempty(self.parent)
                tmp = self.parent.getLabelList();
                labels = {tmp{:} labels{:}};
            end
        end
        
        %%% ------------------------------------------------------
        function plists = getParameterLists(self)
            plists = self.saveItems.parameterLists;
            if ~isempty(self.parent)
                tmp = self.parent.getParameterLists();
                plists = {tmp{:} plists{:}};
            end
        end    
        
        %%% ------------------------------------------------------
        function counts = getParameterCounts(self)
            counts = self.parameterCounts;
            if ~isempty(self.parent)
                tmp = size(self.parent.combinatorics(),1);
                counts = [tmp counts];
            end                  
        end
        
        %%% ------------------------------------------------------
        function lists = getParameterListsFromIdx(self, offset, idx)
            lists = {};
            plists = self.getParameterLists();
            
            for i=1:length(idx)
                if idx(i) <= offset
                    lists{i} = {};
                else
                    lists{i} = plists{idx(i)-offset};
                end
            end
        end
        
        %%% ------------------------------------------------------
        function [idx offset] = parameterNamesToIndices(self, names, varname, wrapfunc)
            idx = zeros(1,length(names));
            labels = self.getLabelList();
            offset = self.getResultVarOffset(varname, wrapfunc);
            for i=1:length(names)
                if isnumeric(names{i})
                    idx(i) = names{i};
                else
                    for k=1:length(labels)
                        if strcmp(labels{k}, names{i})
                            idx(i) = k+offset;
                        end
                    end                    
                    if ~idx(i) 
                        error('Could not find the parameter: %s!', names{i})
                    end
                end
            end            
        end
        
        %%% ------------------------------------------------------
        function idx = getParameterValueIdx(self, offset, nameidx, value)
            constraintParameterList = self.getParameterListsFromIdx(offset, nameidx);
            constraintParameterList=constraintParameterList{1};
            for i=1:length(constraintParameterList)
                if iscell(constraintParameterList) && ischar(constraintParameterList{1})
                    if strcmp(constraintParameterList{i}, value)
                        idx = i;
                        return
                    end
                elseif iscell(constraintParameterList)
                    if constraintParameterList{i} == value
                        idx = i;
                        return
                    end
                elseif constraintParameterList(i) == value
                    idx = i;
                    return                    
                end
            end        
            disp('Value List:')
            disp(constraintParameterList)
            disp('Value:')
            disp(value)
            error('Could not find this value for this parameter!'); 
        end
        
        %%% ------------------------------------------------------
%         function para_vals = getParameterVals(self, idx)
%             if idx > length(self.saveItems.parameterNames)
%                 
%             else
%                 para_vals = self.saveItems.parameterLists{idx};
%             end
%         end
        
        %%% ------------------------------------------------------
        %%% -----------------Plot   FUNCTIONS---------------------
        %%% ------------------------------------------------------
        function plot3(self, varname, parameterNames)
            assert(length(parameterNames)==3, 'You must provide 3 indices, either by their name or index in "varname"!');
            
            M = self.getResultMatrix(varname);
            idx = self.parameterNamesToIndices(parameterNames, varname);
            
            hold on
            classes = unique(M(:,idx(1)));
            for i=1:length(classes)
                subM = M(M(:,idx(1))==classes(i),idx([2 3]));
                plot3(repmat(i,size(subM,1),1),subM(:,1), subM(:,2), '-');
            end 
            labelfuns = {@xlabel, @ylabel, @zlabel};
            for i=1:length(labelfuns)
                if ischar(parameterNames{i})
                    labelfuns{i}(parameterNames{i})
                else
                    labelfuns{i}(sprintf('Idx %d', parameterNames{i}));
                end       
            end
        end
        
        %%% ------------------------------------------------------
        function [ux, uy, Z] = getPlot3Data(self, varname, parameterNames, varargin)
            P.figure = [];
            P.wrapFunc = @(x) x;
            P.constraints = {};
            P.constraintFunc = @(x) x;
            P.constraintName = [];
            P = mysort.util.parseInputs(P, 'getPlot3Data', varargin);            
            M = self.getResultMatrix(varname, P.wrapFunc, P.constraints);
            [idx offset] = self.parameterNamesToIndices(parameterNames, varname, P.wrapFunc);
            
            M = P.constraintFunc(M);
            
            ux = unique(M(:,idx(1)));
            uy = unique(M(:,idx(2)));
            Z = zeros(length(ux), length(uy));
            for i=1:size(M,1)
                Z(ux==M(i,idx(1)), uy==M(i,idx(2))) = M(i,idx(3));
            end            
        end
        
        %%% ------------------------------------------------------
        function [ux, uy, Z] = surf(self, varname, parameterNames, varargin)
            P.figure = [];
            P.wrapFunc = @(x) x;
            P.constraints = {};
            P.constraintFunc = @(x) x;
            P.constraintName = [];
            P = mysort.util.parseInputs(P, 'surf', varargin);
            assert(length(parameterNames)==3, 'You must provide 3 indices, either by their name or index in "varname"!');
            
            [idx offset] = self.parameterNamesToIndices(parameterNames, varname, P.wrapFunc);
            paraLists = self.getParameterListsFromIdx(offset, idx);
            [ux, uy, Z] = self.getPlot3Data(varname, parameterNames, varargin{:});
            
            if isempty(P.figure)
                P.figure = mysort.plot.figure();
            end
            
            surf(ux, uy, Z');
            
            
            labelfuns = {@xlabel, @ylabel, @zlabel};
            ticknames = {'xtick', 'ytick', 'ztick'};
            for i=1:length(labelfuns)
                if ischar(parameterNames{i})
                    labelfuns{i}(parameterNames{i})
                    set(gca, ticknames{i}, 1:length(paraLists{i}));
                    set(gca, [ticknames{i} 'label'], paraLists{i});                    
                else
                    labelfuns{i}(varname);
                end      
            end
            
            if ~isempty(P.constraints)
                str = [];
                for i=1:size(P.constraints,1)
                    str = [str P.constraints{i,1} '=' mysort.util.var2str(P.constraints{i,2})];
                end
                mysort.plot.figureTitle(str);
            elseif ~isempty(P.constraintName)
                mysort.plot.figureTitle(P.constraintName);
            end
        end        
    end
end