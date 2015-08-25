classdef PreProcessedMultiSessionInterface < mysort.ds.ExtendedMultiSessionInterface & ...
                                             mysort.ds.PreProcessedDataSourceInterface
            
    properties
     
    end
    
    methods(Abstract)

    end    
   
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = PreProcessedMultiSessionInterface(preprocessedBufferFilename, preprocessedBufferH5path, h5info, spikeSortingFolder, name, s_per_sec, sessionList)
            assert(isa(sessionList, 'mysort.ds.PreProcessedDataSourceInterface'), 'sessionList must contain objects derived from type mysort.ds.PreProcessedDataSourceInterface!');
            if isempty(name)
                name = 'PreProcessedMultiSessionInterface';
            end
            self = self@mysort.ds.ExtendedMultiSessionInterface(name, s_per_sec,  sessionList);
            self = self@mysort.ds.PreProcessedDataSourceInterface(preprocessedBufferFilename, preprocessedBufferH5path,  h5info, spikeSortingFolder, name, s_per_sec, sessionList(1).getMultiElectrode());
        end
        %------------------------------------------------------------------
        function preprocess(self, sessionIdx)
            if nargin == 1
                sessionIdx = 1:self.getNSessions();
            end
            for i=1:length(sessionIdx)
                self.sessionList(i).preprocess();
            end
        end   
        %------------------------------------------------------------------
        function b = isPreprocessed(self, sessionIdx)
            if nargin == 1
                sessionIdx = self.activeSessionIdx;
            end
            b = self.sessionList(sessionIdx).isPreprocessed();
        end
        %------------------------------------------------------------------            
        function loadSortings(self)
            if isempty(self.spikeSortingFolder)
                self.SpikeSortingContainers = {};
                return
            end
            flist = dir(fullfile(self.spikeSortingFolder, '*.sorting.h5'));
            for i=1:length(flist)
                ffile = fullfile(self.spikeSortingFolder, flist(i).name);
                self.SpikeSortingContainers{end+1} = mysort.mea.PersistentMultiSessionSpikeSortingContainer(ffile, '/SortingContainer');
            end
        end
        %------------------------------------------------------------------
        function smad = noiseStd(self, sessionIndex)
            if nargin < 2
                sessionIndex = self.activeSessionIdx;
            end    
            assert(length(sessionIndex)==1, 'Only one session can be queried at a time!');
            smad = self.sessionList(sessionIndex).noiseStd();
        end

        %------------------------------------------------------------------
        function [times pks] = detectSpikes(self, sessionIndex)
            if nargin < 2
                sessionIndex = self.activeSessionIdx;
            end    
            assert(length(sessionIndex)==1, 'Only one session can be queried at a time!');
            [times pks] = self.sessionList(sessionIndex).detectSpikes();
        end
        %------------------------------------------------------------------
        function smad = getNoiseStdFromBufferFile(self)
            error('Dont call this function');
        end        
        %------------------------------------------------------------------
        function [times pks] = getDetectSpikesFromBufferFile(self)
            error('Dont call this function');
        end  
        %------------------------------------------------------------------
        function MSSSC = createNewSpikeSortingContainer(self, name, gdfList)
            if nargin == 1 || isempty(name)
                name = 'no_name'; 
            end
            fname = fullfile(self.spikeSortingFolder, [name '.sorting.h5']);
            i = 1;
            while exist(fname, 'file')
                fname = fullfile(self.spikeSortingFolder, [name num2str(i) '.sorting.h5']);
                i = i+1;
            end
            if nargin < 3
                gdfList = cell(1, self.getNSessions());
            end
            
            MSSSC = mysort.mea.PersistentMultiSessionSpikeSortingContainer(fname, ...
                name, gdfList);
            MSSSC.save();
            self.SpikeSortingContainers{end+1} = MSSSC;
        end        
    end
end