classdef GridSortInspector < handle
    properties
        %% DEFS
        expName
        path
        configs_path
        files
        h5FileLists
        
        %% LOADED DATA
        configs
        resultsM
        sortings
        
        %% TEMP OBJECTS
        CMOSMEA
        WF
        WFP
        
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = GridSortInspector(expname)
            pd = pdefs();
            self.expName = expname;
            self.configs_path = fullfile(pd.networkTempShare, 'Hillier_2013');
            self.path = fullfile(self.configs_path, [self.expName 'Out']);
            self.files.resultsMichele = fullfile(self.configs_path, [self.expName '_resultsForMichele.mat']);
            self.files.config = fullfile(self.configs_path, [self.expName '.mat']);
            self.configs  = load(self.files.config);
            self.h5FileLists = mysort.mea.configList2H5List(self.configs_path, self.expName, self.configs.configNTKLists);
            
%             self.resultsM = load(self.files.resultsMichele);
            
            self.CMOSMEA = mysort.mea.CMOSMEA(self.h5FileLists{1});
            self.CMOSMEA.concatenateSessions();
            self.loadSortings();
            gdf = self.sortings{1}.gdf_merged;
            gdf = gdf(1:200,:);
%             self.WF = mysort.wf.WfManagerWithDataSource(self.CMOSMEA, gdf(:,2), ones(size(gdf,1),1), gdf(:,1));
%             self.WFP = mysort.plot.WaveformManagerPlot(self.WF, 'plotMedianColor', [.9 .7 .0], 'plotControls', true);
        end
        
        %------------------------------------------------------------------
        function loadSortings(self)
            nConfigs = length(self.h5FileLists);
            for i=1:nConfigs
                sortOutPath = fullfile(self.path, ['Config' num2str(i)], 'Sorting1');
                sourceFile = fullfile(sortOutPath, [self.expName '_proj_' num2str(i) '_results.mat'])
                try
                    assert(exist(sourceFile, 'file')>0, 'Source File not found!');
                    self.sortings{i} = load(sourceFile);
                catch
                end

                % get individual files' length in samples
%                 compoundMea = mysort.mea.compoundMea(self.h5FileLists{i}, 'useFilter', 0, 'name', 'PREFILT');
%                 L = compoundMea.X.getAllSessionsLength();

                % break gdf apart
%                 start_end_times = [0 cumsum(L)];
%                 assert(max(D.gdf_merged(:,2)) <= start_end_times(end), 'Spike Times out of Range !!');
%                 mgdf = mysort.spiketrain.gdf2multiSessionGdf(D.gdf_merged, start_end_times);
%                 for k=1:length(L)
%                     gdf = mgdf(mgdf(:,3)==k,[1:2]);
%                     R{k, i} = gdf;
%                 end
            end            
        end
        
    end
end