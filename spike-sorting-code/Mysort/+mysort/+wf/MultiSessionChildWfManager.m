classdef MultiSessionChildWfManager < mysort.wf.MultiSessionBufferedWfManager
    % This class uses a MultiSessionBufferedWfManager as its WfSource

    properties       
        wfbufferParent
        wfbufferParentEventIdx
        wfbufferParentSessionIdx
    end
    
    methods
        %------------------------------------------------------------------
        function self = MultiSessionChildWfManager(wfbufferParent, subSelectionIdx)       
            self = self@mysort.wf.MultiSessionBufferedWfManager(wfbufferParent.wfDataSource,...
                wfbufferParent.eventTimes(subSelectionIdx), ...
                wfbufferParent.eventChans(subSelectionIdx),...
                wfbufferParent.eventSessions(subSelectionIdx),...
                wfbufferParent.wfbufferCutLeft,...
                wfbufferParent.wfbufferCutLength,...
                wfbufferParent.eventIDs(subSelectionIdx), wfbufferParent.MultiElectrode);
            self.wfbufferParent = wfbufferParent; 
            self.wfbufferParentEventIdx = subSelectionIdx;            
            
            ME = wfbufferParent.MultiElectrode.getSubSessionIdxMultiElectrode(self.wfSessionIdx);            
            [a sidx c] = intersect(wfbufferParent.wfSessionIdx, self.wfSessionIdx);
            self.wfbufferParentSessionIdx = sidx;
        end
        
        %------------------------------------------------------------------
        function setCutLeftAndLength(self, cutleft, cutlength)
            error('you cant do that on a childwfmanager!');
        end
        %------------------------------------------------------------------
        function [wfs subSIdx]= getWaveform4Idx(self, idx, varargin)
            parentIdx = self.wfbufferParentEventIdx(idx);
            wfs = self.wfbufferParent.getWaveform4Idx(parentIdx, varargin{:});
            subSIdx = self.wfbufferParentSessionIdx;
            wfs = wfs(subSIdx);
        end
        %------------------------------------------------------------------
        function wfManager = getWfManger4SubIdx(self, idx)
            wfManager = mysort.wf.MultiSessionChildWfManager(...
                self.wfbufferParent, self.wfbufferParentEventIdx(idx));
        end          
    end
end