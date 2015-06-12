classdef SpikeSortingContainer < handle
    properties
        name
        gdf
        wfDataSource
        templateWfs % matrix, time x channels x units
        templateCutLeft
        templateCutLength
        unitNames
        tMinMax
        nMaxSpikesForTemplateCalc
        bIsGroundTruth
    end
    
    
    methods
        %------------------------------------------------------------------
        function self = SpikeSortingContainer(name, gdf, varargin)
            P.wfDataSource = [];
            P.nMaxSpikesForTemplateCalc = 100;
            if isstruct(name)
                self.fromStruct(name);
                P = mysort.util.parseInputs(P, varargin, 'error');
                self.wfDataSource = P.wfDataSource;
                return
            end
            
            P.unitNames = [];
            P.isGroundTruth = false;
            P.templateWfs = [];
            P.templateCutLeft = [];
            P.templateCutLength = [];
            P = mysort.util.parseInputs(P, varargin, 'error');
            self.name = name;
            if iscell(gdf)
                self.gdf = mysort.spiketrain.toGdf(gdf);
            else
                self.gdf = gdf;
            end
            if isempty(P.unitNames) && ~isempty(gdf)
                P.unitNames = unique(gdf(:,1));
            end
            self.unitNames = P.unitNames;
            if ~isempty(gdf)
                self.tMinMax = [min(gdf(:,2)) max(gdf(:,2))];
            else
                self.tMinMax = [0 1];
            end
            self.bIsGroundTruth = P.isGroundTruth;
            self.templateWfs = P.templateWfs;
            if ~isempty(P.templateWfs)
                assert(~isempty(P.templateCutLeft), 'If template waveforms are given, cutleft has to be set!');
                if ~isempty(P.templateCutLength)
                    assert(P.templateCutLength == size(self.templateWfs,1), 'TemplateCutLength does not match given waveforms!');
                end
            end
            self.templateCutLeft = P.templateCutLeft;
            self.templateCutLength = P.templateCutLength;
            self.wfDataSource = P.wfDataSource;
            self.nMaxSpikesForTemplateCalc = P.nMaxSpikesForTemplateCalc;
        end
        %------------------------------------------------------------------
        function fromStruct(self, S)
            self.name = S.name;
            self.gdf = S.gdf;
            self.unitNames = S.details.unitNames;
            self.tMinMax = S.details.tMinMax;
            self.bIsGroundTruth = S.details.isGroundTruth;
            if isempty(S.details.templateWfs)
                self.templateWfs = [];
            else
                nC = S.details.templateNC;
                self.templateWfs = mysort.wf.v2t(S.details.templateWfs, nC);
            end
            self.templateCutLeft = S.details.templateCutLeft;            
        end
        %------------------------------------------------------------------
        function S = toStruct(self)
            S.name = self.name;
            S.gdf = self.gdf;
            S.details.unitNames = self.unitNames;
            S.details.tMinMax = self.tMinMax;
            S.details.isGroundTruth = self.bIsGroundTruth;
            S.details.templateNC = size(self.templateWfs, 2);
            S.details.templateWfs = mysort.wf.t2v(self.templateWfs);
            S.details.templateCutLeft = self.templateCutLeft;
            S.date_ = date();
            S.version = 1;
            S.readme = 'This is a struct derived from the class SpikeSortingContainer. Dont edit if you dont know what you are doing.';
        end        
        %------------------------------------------------------------------
        function save2File(self, fname, h5path)
            S = self.toStruct();
            mysort.h5.recursiveSave(fname, S, h5path);
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
        function idx = unitNames2UnitIdx(self, unitNames, varargin)
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
        function gdf = getGdf(self, t1, t2, unitNames)
            if nargin < 4 || isempty(unitNames)
                unitIdx = 1:length(self.unitNames);
            else
                unitIdx = self.unitNames2UnitIdx(unitNames);
            end
            if nargin == 1
                gdf = self.getGdf4UnitIdx(unitIdx);
            else
                gdf = self.getGdf4UnitIdx(unitIdx, t1, t2);
            end
        end
        
        %------------------------------------------------------------------
        function gdf = getGdf4UnitIdx(self, unitIdx, t1, t2)
            if isempty(self.gdf)
                gdf = [];
                return
            end
            if nargin == 1
                gdf = self.gdf;
                return
            end
            if nargin > 2
                gdf = self.gdf(self.gdf(:,2) <= t2 & self.gdf(:,2) >= t1,:); 
            else
                gdf = self.gdf;
            end
            if ~isempty(unitIdx)
                gdf = gdf(ismember(gdf(:,1), self.unitNames(unitIdx)),:);
            end
        end
        %------------------------------------------------------------------
        function T = getTemplateWaveforms(self, varargin)
            unitIdx = self.unitNames2UnitIdx(varargin{:});
            T = self.getTemplateWaveforms4UnitIdx(unitIdx);
        end
        %------------------------------------------------------------------
        function T = getTemplateWaveforms4UnitIdx(self, unitIdx)
            if isempty(self.templateWfs)
                self.computeTemplateWfs();
            end
            T = self.templateWfs(:, :, unitIdx);
        end
        %------------------------------------------------------------------
        function T = computeTemplateWfs(self, cutLeft, cutLength)
            assert(~isempty(self.wfDataSource), 'tempates are requested, but were not given, and no datasource is available to compute them!');
            if nargin < 2
                cutLeft = self.templateCutLeft;
            end
            if isempty(cutLeft)
                cutLeft = 20;
            end
            
            if nargin < 3
                cutLength = self.templateCutLength;
            end
            if isempty(cutLength)
                cutLength = 80;
            end
            nC = self.wfDataSource.MultiElectrode.getNElectrodes();
            T = zeros(cutLength, nC, length(self.unitNames));
            
            disp('computing templates, that might take a while...')
            for i=1:length(self.unitNames)
                fprintf('%d of %d\n', i, length(self.unitNames));
                wfs = self.getWaveforms4UnitIdx(i, cutLeft, cutLength, self.nMaxSpikesForTemplateCalc);
                if isempty(wfs)
                    continue
                end
                T(:,:,i) = mysort.wf.v2m(median(wfs,1), nC)';
            end
            self.templateCutLeft = cutLeft;
            self.templateCutLength = cutLength;
            self.templateWfs = T;
        end
        
        %------------------------------------------------------------------
        function TL = getTemplateList(self)
            WFS = self.getTemplateWaveforms();
            TL = mysort.mea.Template.empty();
            for n = 1:size(WFS,3);
                TL(n) = mysort.mea.Template(WFS(:,:,n), self.templateCutLeft, ...
                    self.wfDataSource.getSamplesPerSecond(), n, self.wfDataSource.MultiElectrode.copy(), ones(1, self.wfDataSource.MultiElectrode.getNElectrodes()));
            end
        end
            
        %------------------------------------------------------------------
        function wfs = getWaveforms4UnitIdx(self, i, cutLeft, cutLength, Nmax)
            if nargin < 5
                Nmax = 100;
            end
            t = self.getGdf4UnitIdx(i);
            wfs = [];
            if isempty(t)
                return
            end
            if size(t,1) > Nmax
                t = t(1:Nmax,2);
            else
                t = t(:,2);
            end
            % for the moment round t!
            t = round(t);
            wfs = self.wfDataSource.getWaveform(t, cutLeft, cutLength);
        end
        
        %------------------------------------------------------------------
        function cl = getTemplateCutLeft(self)
            cl = self.templateCutLeft;
        end
        %------------------------------------------------------------------
        function SSC = getSubSpikeSortingContainer(self, unitNames, varargin)
            if isempty(unitNames)
                unitIdx = 1:length(self.unitNames);
            else
                unitIdx = self.unitNames2UnitIdx(unitNames);
            end
            SSC = self.getSubSpikeSortingContainer4UnitIdx(unitIdx, varargin{:});
        end
        %------------------------------------------------------------------
        function SSC = getSubSpikeSortingContainer4UnitIdx(self, unitIdx, varargin)
            gdf = self.getGdf4UnitIdx(unitIdx, varargin{:});
            unitNames = self.unitNames(unitIdx);
            SSC = mysort.spiketrain.SpikeSortingContainer(...
                self.name, gdf, 'unitNames', unitNames,...
                'templateWfs', self.templateWfs, ...
                'templateCutLeft', self.templateCutLeft,...
                'isGroundTruth', self.bIsGroundTruth);
        end       
    end
    methods(Static)
        %------------------------------------------------------------------
        function S = emptyStruct()
            S.name = '';
            S.gdf = [];
            S.unitNames = [];
            S.tMinMax = [];
            S.isGroundTruth = false;
            S.templateWfs = [];
            S.templateCutLeft = 0;
            S.date_ = date();
            S.version = 1;
            S.readme = 'This is a struct derived from the class SpikeSortingContainer. Dont edit if you dont know what you are doing.';
        end       
    end
end