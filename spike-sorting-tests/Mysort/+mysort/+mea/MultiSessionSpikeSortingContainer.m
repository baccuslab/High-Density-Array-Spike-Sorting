classdef MultiSessionSpikeSortingContainer < handle
    properties
        name
        
        singleSessionGdfList
        TemplateManager
        unitNames
        bIsGroundTruth   
    end
    
    
    methods
        %------------------------------------------------------------------
        function self = MultiSessionSpikeSortingContainer(name, gdfList, varargin)
            % This function can be called in three different ways:
            % 1. SC = mysort.mea.MultiSessionSpikeSortingContainer(S)
            %    where S is a structure ceared by to "toStruct" function
            %    recreates a e.g. stored container
            % 2. SC = mysort.mea.MultiSessionSpikeSortingContainer(name, TM,
            %            cutleft, singleSessionGdfList
            %    NEEDS TO BE IMPLEMENTED
            % 3. SC = mysort.mea.MultiSessionSpikeSortingContainer(name, gdfList)
            %    Computes all necessary information itself.
            
            if isstruct(name)
                self.fromStruct(name);
                return
            end
            P.unitNames = [];
            P.isGroundTruth = false;
            P = mysort.util.parseInputs(P, varargin, 'error');
            if isnumeric(gdfList)
                % is a gdf, but multisession?
                if size(gdfList,2) == 2
                    % yes, classical 2 column gdf
                    gdfList = {gdfList};
                elseif size(gdfList,2) == 3
                    % no, this is a multisession gdf, with a third column
                    sessions = unique(gdfList(:,3));
                    gdf = {};
                    for sidx=1:length(sessions)
                        s = sessions(sidx);
                        gdf{s} = gdfList(gdfList(:,3)==s,1:2);
                    end
                    gdfList = gdf;
                else
                    error('unkown second parameter!');
                end
            end
            self.singleSessionGdfList = gdfList;

            if isempty(P.unitNames)
                for i=1:length(gdfList)
                    gdf = gdfList{i};
                    if ~isempty(gdf)
                        P.unitNames = unique([P.unitNames(:); unique(gdf(:,1))]);
                    end
                end
            end
            self.name = name;
            self.unitNames = P.unitNames;
            self.bIsGroundTruth = P.isGroundTruth;
        end
        %------------------------------------------------------------------
        function fromStruct(self, S)
            self.name = S.name;
            if isfield(S, 'singleSessionGdfList')
                self.singleSessionGdfList = S.singleSessionGdfList;
            else
                self.singleSessionGdfList = {};
            end
            self.unitNames = S.details.unitNames;
            self.bIsGroundTruth = S.details.isGroundTruth;
            if isempty(S.TemplateManager)
                self.TemplateManager = [];
            else
                self.TemplateManager = mysort.mea.TemplateManager(S.TemplateManager);
            end
        end
        %------------------------------------------------------------------
        function S = toStruct(self)
            S.name = self.name;
            S.singleSessionGdfList = self.singleSessionGdfList;
            S.details.unitNames = self.unitNames;
            S.details.isGroundTruth = self.bIsGroundTruth;
            if isempty(self.TemplateManager)
                S.TemplateManager = [];
            else
                S.TemplateManager = self.TemplateManager.toStruct();
            end
            S.date_ = date();
            S.version = 1;
            S.readme = 'This is a struct derived from the class MultiSessionSpikeSortingContainer. Dont edit if you dont know what you are doing.';
        end        
        %------------------------------------------------------------------
        function name = getName(self)
            name = self.name;
        end
        %------------------------------------------------------------------
        function b = isGroundTruth(self)
            b = self.bIsGroundTruth;
        end
        
        %------------------------------------------------------------------
        function idx = unitNames2UnitIdx(self, unitNames)
            if nargin == 1 || isempty(unitNames)
                idx = 1:length(self.unitNames);
                return
            end
            idx = zeros(1, length(unitNames));
            for i=1:length(unitNames)
                idx(i) = find(self.unitNames == unitNames(i),1);
            end
        end
        %------------------------------------------------------------------
        function gdf = getGdf(self, sessionidx, t1, t2, unitNames)
            if nargin < 5 || isempty(unitNames)
                unitIdx = 1:length(self.unitNames);
            else
                unitIdx = self.unitNames2UnitIdx(unitNames);
            end
            if nargin == 2
                gdf = self.getGdf4UnitIdx(sessionidx, unitIdx);
            else
                gdf = self.getGdf4UnitIdx(sessionidx, unitIdx, t1, t2);
            end
        end
        
        %------------------------------------------------------------------
        function gdf = getGdf4UnitIdx(self, sessionidx, unitIdx, t1, t2)
            gdf = self.singleSessionGdfList{sessionidx};
            if nargin == 2
                return
            end
            if nargin > 3
                gdf = gdf(gdf(:,2) <= t2 & gdf(:,2) >= t1,:); 
            end
            if nargin > 2 && ~isempty(unitIdx)
                removeUnits = setdiff(1:length(self.unitNames), unitIdx);
                for i=1:length(removeUnits)
                    gdf(gdf(:,1) == removeUnits(i),:) = [];
                end
            end
        end
        %------------------------------------------------------------------
        function gdf = getMultiSessionGdf(self, sessionidx, unitIdx, t1, t2)
            gdf = [];
            if isempty(sessionidx)
                return
            end
            for i=1:length(sessionidx)
                if sessionidx(i) <= length(self.singleSessionGdfList)
                    sgdf = self.singleSessionGdfList{sessionidx(i)};
                    gdf = [gdf ; [sgdf sessionidx(i)*ones(size(sgdf,1),1)]];
                end
            end
            if nargin > 3
                gdf = gdf(gdf(:,2) <= t2 & gdf(:,2) >= t1,:); 
            end
            if nargin > 2 && ~isempty(unitIdx)
                removeUnits = setdiff(1:length(self.unitNames), unitIdx);
                for i=1:length(removeUnits)
                    gdf(gdf(:,1) == removeUnits(i),:) = [];
                end
            end
        end        
        %------------------------------------------------------------------
        function S = getSingleSessionSorting(self, sessionidx, ME)
            if ~isempty(self.TemplateManager)
                [wfs cutleft] = self.TemplateManager.getWaveforms4MultiElectrode(ME);
                S = mysort.spiketrain.SpikeSortingContainer([self.name '_child'], self.getGdf(sessionidx), ...
                    'unitNames', self.unitNames, 'templateWfs', wfs, 'templateCutLeft', cutleft);                
            else
                S = mysort.spiketrain.SpikeSortingContainer([self.name '_child'], ...
                    self.getGdf(sessionidx), 'unitNames', self.unitNames);   
            end
        end
        %------------------------------------------------------------------
        function computeTemplates(self, wfDataSource, varargin)
            P.maxNWaveforms = 1000;
            P = mysort.util.parseInputs(P, varargin);
            
            if isempty(self.TemplateManager) || isempty(self.TemplateManager.TemplateList)
                self.TemplateManager = mysort.mea.TemplateManager(mysort.mea.Template.empty());
                joinedgdf = cell2mat(self.singleSessionGdfList(~cellfun(@isempty, self.singleSessionGdfList))');
                if isempty(joinedgdf)
                    tNames = [];
                    nT = 0;
                else
                    tNames = unique(joinedgdf(:,1));
                    nT = length(tNames);
                end
                cutLeft = 10;
                cutLength = 50;
            else
                cutLeft = self.TemplateManager.TemplateList(1).cutLeft;
                cutLength = 80;
                nT = length(self.TemplateManager.TemplateList);
                tNames = [self.TemplateManager.TemplateList.name];
            end
            for i=1:nT
                ts = [];
                sessions = [];
                for s=1:length(self.singleSessionGdfList)
                    mySgdf = self.singleSessionGdfList{s};
                    if ~isempty(mySgdf)
                        mySgdf = mySgdf(mySgdf(:,1)== tNames(i),2);
                        ts = [ts; mySgdf(:)];
                        sessions = [sessions; s*ones(length(mySgdf),1)];
                    end
                end
                if ~isempty(P.maxNWaveforms) && length(ts) > P.maxNWaveforms 
                    rp = randperm(length(ts));
                    rp = rp(1:P.maxNWaveforms);
                else
                    rp = 1:length(ts);
                end
                    
                WFM = mysort.wf.MultiSessionBufferedWfManager(wfDataSource,...
                    ts, i*ones(length(ts),1), sessions, ...
                    cutLeft, cutLength, 1:length(ts), wfDataSource.getAllSessionsMergedMultiElectrode());
                [T, ME, wfsAllSessions, nWfsPerElectrode] = WFM.getTemplate4Idx(rp);
                if 0
                    figure
                    ax = axes;
                    hold on;
                    for s=1:length(wfsAllSessions)
                        EP = ME.MultiElectrodeList(s).electrodePositions;
                        wfs = mysort.wf.v2t(wfsAllSessions{s}, ME.MultiElectrodeList(s).getNElectrodes());
                        
                        mysort.plot.waveforms2D(wfs, EP, 'plotArgs', {'color', [.7 .7 .8]}, 'AxesHandle', ax);
                    end
                    mysort.plot.waveforms2D(T', ME.electrodePositions, 'plotArgs', {'color', 'r'}, 'AxesHandle', ax);
                end
                if length(self.TemplateManager.TemplateList)>=i
                    self.TemplateManager.TemplateList(i).waveforms = T';
                    self.TemplateManager.TemplateList(i).nSourceSpikesPerElectrode = nWfsPerElectrode;
                    self.TemplateManager.TemplateList(i).MultiElectrode = ME;
                else
                    template = mysort.mea.Template(T', cutLeft, wfDataSource.samplesPerSecond, tNames(i), ME, nWfsPerElectrode);
                    self.TemplateManager.TemplateList(i) = template;
                end
            end
        end
        
        %------------------------------------------------------------------
        function addTemplateList(self, singleSessionGdfList, TemplateList)
            for i=1:length(TemplateList)
                idx = find(self.unitNames == TemplateList(i).name,1);
                assert(isempty(idx), 'Template is already in Template Manager!');
                self.unitNames(end+1) = TemplateList(i).name;
            end
            
            if ~isempty(self.singleSessionGdfList)
                assert(length(singleSessionGdfList) == length(self.singleSessionGdfList), 'the gdf list must have one entry per session!');
                for i=1:length(self.singleSessionGdfList)
                    if length(self.singleSessionGdfList) < i || isempty(self.singleSessionGdfList{i})
                        self.singleSessionGdfList{i} = singleSessionGdfList{i};
                    elseif ~isempty(singleSessionGdfList{i})
                        self.singleSessionGdfList{i} = sortrows(...
                            [self.singleSessionGdfList{i}; singleSessionGdfList{i}], 2);
                    end
                end
                if isempty(self.TemplateManager)
                    self.TemplateManager = mysort.mea.TemplateManager(TemplateList);
                else
                    self.TemplateManager.TemplateList = [self.TemplateManager.TemplateList TemplateList];
                end
            else
                self.singleSessionGdfList = singleSessionGdfList;
                self.TemplateManager = mysort.mea.TemplateManager(TemplateList);
            end            
        end
        %------------------------------------------------------------------
        function deleteTemplate(self, unitNames)
            if isempty(unitNames)
                return
            end
            idx = self.unitNames2UnitIdx(unitNames(:));
            for i=1:length(idx)
                for s = 1:length(self.singleSessionGdfList)
                    S = self.singleSessionGdfList{s};
                    S(S(:,1)==unitNames(i),:) = [];
                    self.singleSessionGdfList{s} = S;
                end
            end
            self.TemplateManager.deleteTemplateIdx(idx);
            self.unitNames(idx) = [];
        end
    end
end