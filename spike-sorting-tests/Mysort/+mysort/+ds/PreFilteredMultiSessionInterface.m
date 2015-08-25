classdef PreFilteredMultiSessionInterface < mysort.ds.FilteredMultiSessionInterface & ...
                                            mysort.ds.PreFilteredDataSourceInterface
            
    properties
     
    end
    
    methods(Abstract)

    end    
   
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = PreFilteredMultiSessionInterface(prefilterBufferFilename, prefilterBufferH5path, h5info, filterFactory, useFilter, name, s_per_sec, sessionList)
            assert(isa(sessionList, 'mysort.ds.PreFilteredDataSourceInterface'), 'sessionList must contain objects derived from type mysort.ds.PreFilteredDataSourceInterface!');
            if isempty(name)
                name = 'PreFilteredMultiSessionInterface';
            end
            self = self@mysort.ds.FilteredMultiSessionInterface(filterFactory, useFilter, name, s_per_sec, sessionList);
            self = self@mysort.ds.PreFilteredDataSourceInterface(prefilterBufferFilename, prefilterBufferH5path, h5info, filterFactory, useFilter, name, s_per_sec, sessionList(1).getMultiElectrode());
        end
        %------------------------------------------------------------------
        function prefilter(self, sessionIdx)
            if nargin == 1
                sessionIdx = 1:self.getNSessions();
            end
            for i=1:length(sessionIdx)
                self.sessionList(i).prefilter();
            end
        end   

        %------------------------------------------------------------------
        function X = getFilterDataFromBufferFile(self, varargin)
            error('Dont call this method!');
        end           
    end
end