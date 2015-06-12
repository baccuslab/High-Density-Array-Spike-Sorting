classdef God < handle
    properties
        DataSource
        SortingContainer
        templateCutLeft
        templateCutLength
        wfBuffer
    end
    
    
    methods   
        %------------------------------------------------------------------
        function self = God(DataSource, SortingContainer, varargin)
            P.templateCutLeft = 10;
            P.templateCutLength = 60;
            P = mysort.util.parseInputs(P, varargin, 'error');
            
            self.DataSource = DataSource;
            self.setSortingContainer(SortingContainer);
            self.wfBuffer = [];
            self.templateCutLeft = P.templateCutLeft;
            self.templateCutLength = P.templateCutLength;
        end
       
        %------------------------------------------------------------------
        function setSortingContainer(self, SC)
            self.SortingContainer = SC;
        end

        %------------------------------------------------------------------
        function T = getTemplateWaveforms(self, sessionidx, unitNames, nMaxUseWaveforms)
            if nargin < 3 || isempty(unitNames)
                unitIdx = self.SortingContainer.unitNames2UnitIdx();
            else
                unitIdx = self.SortingContainer.unitNames2UnitIdx(unitNames);
            end
            if nargin < 4 
                nMaxUseWaveforms = 25;
            end
            T = self.getTemplateWaveforms4UnitIdx(sessionidx, unitIdx, nMaxUseWaveforms);
        end
        %------------------------------------------------------------------
        function T = getTemplateWaveforms4UnitIdx(self, sessionidx, unitIdx, nMaxUseWaveforms)
            if nargin < 4 
                nMaxUseWaveforms = 25;
            end            
            if isempty(self.wfBuffer)
                self.initWfBuffers();
            end
            T = mysort.mea.Template.empty();
            for i=1:length(unitIdx)
                wfidx = 1:min(self.wfBuffer(unitIdx(i)).getNWfs(), nMaxUseWaveforms);
                [t, ME, wfsAllSessions, nWfs] = self.wfBuffer(unitIdx(i)).getTemplate4Idx(wfidx);
                T(i) = mysort.mea.Template(mysort.wf.m2v(t), self.templateCutLeft,...
                    self.DataSource.getSamplesPerSecond, unitIdx(i), ME, nWfs);
            end
        end
        %------------------------------------------------------------------
        function initWfBuffers(self)
            mgdf = mysort.spiketrain.mergeGdfs(self.SortingContainer.singleSessionGdfList);
            self.wfBuffer = mysort.wf.MultiSessionBufferedWfManager.empty();
            ME = self.DataSource.getAllSessionsMergedMultiElectrode();
            for i=1:length(self.SortingContainer.unitNames)
                mgdfidx = mgdf(:,2)==self.SortingContainer.unitNames(i);
                self.wfBuffer(i) = mysort.wf.MultiSessionBufferedWfManager(...
                    self.DataSource,...
                    mgdf(mgdfidx,3), mgdf(mgdfidx,2), mgdf(mgdfidx,1), ...
                    self.templateCutLeft, self.templateCutLength,...
                    1:sum(mgdfidx), ME);
            end
        end
        %------------------------------------------------------------------
        function wfs = getWaveforms4UnitIdx(self, i, cutLeft, cutLength, Nmax)

        end
        
        %------------------------------------------------------------------
        function cl = getTemplateCutLeft(self)
            cl = self.templateCutLeft;
        end    
    end
end