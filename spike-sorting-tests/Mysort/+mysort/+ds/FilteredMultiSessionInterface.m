classdef FilteredMultiSessionInterface < mysort.ds.FilteredDataSourceInterface &...
                                         mysort.ds.MultiSessionInterface        
            
    properties
     
    end
    
    methods(Abstract)

    end    
   
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = FilteredMultiSessionInterface(filterFactory, useFilter, name, s_per_sec, sessionList)
            assert(isa(sessionList, 'mysort.ds.PreFilteredDataSourceInterface'), 'sessionList must contain objects derived from type mysort.ds.PreFilteredDataSourceInterface!');
            if isempty(name)
                name = 'PrefilteredMultiSessionInterface';
            end
            self = self@mysort.ds.FilteredDataSourceInterface(filterFactory, useFilter, name, s_per_sec, sessionList(1).getMultiElectrode());
            self = self@mysort.ds.MultiSessionInterface(name, s_per_sec, sessionList);            
        end
       %------------------------------------------------------------------
        function X = getUnfilteredData(self, timeIndex, channelIndex, sessionIndex, varargin)
            if nargin < 4
                sessionIndex = self.activeSessionIdx;
                if nargin < 3
                    channelIndex = 1:self.size_(2);
                    if nargin < 2
                        timeIndex = 1:self.size_(1);
                    end
                end
            end    
            X = self.sessionList(sessionIndex).getUnfilteredData(timeIndex, channelIndex, varargin{:});
        end       
    end
    methods(Sealed)
        %------------------------------------------------------------------
        function X = getFilteredData(self, timeIndex, channelIndex, sessionIndex, varargin)
            if nargin < 4
                sessionIndex = self.activeSessionIdx;
                if nargin < 3
                    channelIndex = 1:self.size_(2);
                    if nargin < 2
                        timeIndex = 1:self.size_(1);
                    end
                end
            end    
            X = self.sessionList(sessionIndex).getFilteredData(timeIndex, channelIndex, varargin{:});
        end      
    end
end
 