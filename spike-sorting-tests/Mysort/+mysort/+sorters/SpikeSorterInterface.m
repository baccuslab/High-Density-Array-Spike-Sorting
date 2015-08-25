
classdef SpikeSorterInterface  < mysort.util.DebuggableClass
    properties
        bReadyToSort
        bSorted 
          
        DH
        Tf 
        
        sorting
        templates 
    end
    methods (Abstract)
        sorting = sort_(self)
    end    
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = SpikeSorterInterface()  
            self = self@mysort.util.DebuggableClass();
           	self.DH = [];
            self.Tf = [];
            self.bReadyToSort = false;        
            self.bSorted = false;

            self.sorting = [];
            self.templates = [];
        end   

        %%% ------------------------------------------------------ 
        function sorting = sort(self, X)
            if ~isa(X, 'mysort.ds.DataSourceInterface')
                self.checkX(X);
                X = mysort.ds.Matrix(X');
            end
            self.checkDataFile(X);
            self.DH = X;

            sorting = self.sort_();
            self.sorting = sorting;
            self.bSorted = true;
%             if self.P.cacheable
%                 self.cache_me();
%             end
        end
        %%% ------------------------------------------------------ 
        function checkX(self, X)
            assert(ndims(X)==2, 'Data needs to be a matrix containing the recording channels as rows');
            assert(size(X,2)>size(X,1), 'Data needs to be transposed? More channels than samples...'); 
        end
        %%% ------------------------------------------------------ 
        function checkDataFile(self, df)
            assert(~isempty(df), 'If sort is called without X, you need to provide a DataSourceInterface instance!');
            assert(isa(df, 'mysort.ds.DataSourceInterface'), 'Not a DataSourceInterface instance!');
        end

      
        %%% ------------------------------------------------------
        function [artefactEpochs signal threshold] = computeArtefactEpochs(self, X)
%             [artefactEpochs signal threshold] =  zeroCrossingArtefactDetector(X);
            artefactEpochs = [];
            signal = [];
            threshold = [];
        end
        
        %%% ------------------------------------------------------
        function [X start stopp] = getX(self, X, start, stopp)
            error('This function is obsolete! Use the one from self.DH !');
        end

        %%% ------------------------------------------------------
        function templates = getTemplates(self)
            templates = self.templates;
        end               
        %%% ------------------------------------------------------
        function classes = getSpikeClasses(self, varargin)
            P.start = [];
            P.stopp = [];
            P.restrict = [];
            P = mysort.util.parseInputs(P, 'getSpikeClasses', varargin);  
            if isempty(P.start); P.start = 1; end
            if isempty(P.stopp); P.stopp = size(self.DH,1); end
            classes = self.sorting(self.sorting(:,2)>=P.start & self.sorting(:,2)<=P.stopp,1);
            if ~isempty(P.restrict)
                classes = classes(ismember(classes, P.restrict));
            end
        end        
        %%% ------------------------------------------------------
        function spikeEpochs = getSpikeEpochs(self, varargin)
            P.start = [];
            P.stopp = [];
            P.restrict = [];
            P = mysort.util.parseInputs(P, 'getSpikeEpochs', varargin);
            if isempty(P.start); P.start=1; end
            if isempty(P.stopp); P.stopp=self.DH.Len; end
            if isempty(self.sorting)
                spikeEpochs = [];
                return
            end
            idx = self.sorting(:,2)>=P.start & self.sorting(:,2)<=P.stopp;
            if ~isempty(P.restrict)
                idx(idx) = idx(idx) & ismember(self.sorting(idx,1), P.restrict);
            end
            spikeEpochs = [self.sorting(idx,2) self.sorting(idx,2)+self.Tf-1];
        end                
        %%% ------------------------------------------------------
        function [spikes, classes] = getSpikeWaveforms(self, varargin)
            P.restrict = [];
            P.start = [];
            P.stopp = [];
            P = mysort.util.parseInputs(P, 'getSpikeWaveforms', varargin);    
            if isempty(P.start); P.start = 1; end;
            classes = self.getSpikeClasses('start', P.start,'stopp', P.stopp, 'restrict', P.restrict);
            epochs = self.getSpikeEpochs('start', P.start, 'stopp', P.stopp, 'restrict', P.restrict);
            X = self.DH(P.start:P.stopp,:)';
            spikes = mysort.epoch.extractWaveform(X, epochs - P.start +1);
        end
        
        %%% ------------------------------------------------------
        function gdf = getSorting(self,varargin)
            P.restrict = [];
            P.start = [];
            P.stopp = []; 
            P = mysort.util.parseInputs(P, 'getSpikeWaveforms', varargin); 
            if isempty(P.start); P.start = 1; end
            if isempty(P.stopp); P.stopp = size(self.DH,1); end         
            if isempty(self.sorting)
                gdf = [];
                return
            end
            idx = self.sorting(:,2)>=P.start & self.sorting(:,2)<=P.stopp;
            if ~isempty(P.restrict)
                idx(idx) = idx(idx) & ismember(self.sorting(idx,1), P.restrict);
            end            
            gdf = self.sorting(idx,:);
        end
        
        %%% ------------------------------------------------------
        function gtgdf = getGtGdf(self, R, varargin)
            P.restrict = [];
            P.start = [];
            P.stopp = []; 
            P = mysort.util.parseInputs(P, 'getGTSpikeTrain', varargin); 
            if isempty(P.start); P.start = 1; end
            if isempty(P.stopp); P.stopp = size(self.DH,1); end  
            gtgdf = mysort.spiketrain.toGdf(R.St1);
            idx = gtgdf(:,2)>=P.start & gtgdf(:,2)<=P.stopp;
            if ~isempty(P.restrict)
                idx(idx) = idx(idx) & ismember(gtgdf(idx,1), P.restrict);
            end            
            gtgdf = gtgdf(idx,:);
        end
        
        %%% ------------------------------------------------------
        function plotSpikes(self, varargin)
            [spikes classes] = self.getSpikeWaveforms(varargin{:});
            mysort.plot.spikes(spikes, 'classes', classes, 'nC', size(elf.DH,2));
            mysort.plot.figureName('Spikes');                        
        end         
        %%% ------------------------------------------------------
        function plotTemplates(self)
            mysort.plot.spikes(self.templates, ...
                        'classes', 1:size(self.templates,1),...
                        'nC', size(self.DH,2), 'linewidth', 3);
            mysort.plot.figureName('Templates');                        
        end   
        %%% ------------------------------------------------------
        function plotClusterPCA(self, varargin)
            [spikes classes] = self.getSpikeWaveforms(varargin{:});
            fetX = mysort.util.dimReductionPCA(spikes, 4);
            mysort.plot.clustering(fetX, classes);
            mysort.plot.figureName('ClusterPCA'); 
        end

        %%% ------------------------------------------------------
        function plotISI(self,varargin)
            gdf = self.getSorting(varargin{:});
            mysort.plot.isi(gdf, 'print2ms', 1, 'isHist', 0, 'srate', self.DH.getSamplesPerSecond())
            mysort.plot.figureName('ISI'); 
        end
        
        %%% ------------------------------------------------------
        function P = plotSorting(self, varargin)
            import mysort.*
            P.start = [];
            P.stopp = [];
            P.axesHandle = [];
            P.figHandle = [];
            P.spacer = [];
            P.restrict = [];
            P.gtSpikeTrain = [];
            P.eval = [];
            P.cutleft = 0;
            P = util.parseInputs(P, 'plotSorting', varargin, 1);
            if ~isempty(P.stopp); P.stopp = min(P.stopp, size(self.DH,1)); end
            X = self.DH(P.start:P.stopp,:)';
            gdf = self.getSorting('start', P.start, 'stopp', P.stopp, 'restrict', P.restrict);
            gtgdf = [];
            if ~isempty(P.eval)
                gtgdf = self.getGtGdf(P.eval, 'start', P.start, 'stopp', P.stopp, 'restrict', P.restrict);
            end
            mysort.plot.sorting(X, gdf, mysort.wf.v2t(self.templates, size(X,1)), 'cutleft', P.cutleft, ...
                'spacer', P.spacer, 'sampleOffset', P.start, 'gtGdf', gtgdf, 'template_IDs', 1:size(self.templates,1))
        end
    end
end