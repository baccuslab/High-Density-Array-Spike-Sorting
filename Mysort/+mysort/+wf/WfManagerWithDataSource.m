classdef WfManagerWithDataSource < mysort.wf.WfManager
    properties             
        wfDataSource
    end
    
    %------------------------------------------------------------------
    methods (Abstract)

    end
    
    methods
        %------------------------------------------------------------------
        function self = WfManagerWithDataSource(wfDataSource, eventTimes, eventChans, eventIDs, ME)
            if nargin < 5 || isempty(ME)
                ME = wfDataSource.MultiElectrode;
            end
            self = self@mysort.wf.WfManager(eventTimes, eventChans, eventIDs, ME);
            self.wfDataSource = wfDataSource;

        end
        %------------------------------------------------------------------
        function wfs = getWaveform4Idx(self, eventIdx, cutleft, cutlength, electrodeNumber)
            if nargin < 5
                electrodeNumber = [];
            end            
            if nargin < 3
                cutleft = 10;
            end
            if nargin < 4
                cutlength = 60;
            end
            elIdx = self.MultiElectrode.getElIdx4ElNumber(electrodeNumber);
            t = self.eventTimes(eventIdx);
            if 0
                %%
                id = self.eventChans(eventIdx);
            end
            wfs = self.wfDataSource.getWaveform(t, cutleft, cutlength, elIdx);
        end
        
        %------------------------------------------------------------------
        function wfManager = getWfManger4SubIdx(self, idx)
            wfManager = mysort.wf.WfManager(self.wfDataSource,...
                self.eventTimes(idx), self.eventChans(idx), self.eventIDs(idx));
        end
        
        %------------------------------------------------------------------
        function nC = getNChannels(self)
            nC = self.wfDataSource.MultiElectrode.getNElectrodes();
        end
    end
end