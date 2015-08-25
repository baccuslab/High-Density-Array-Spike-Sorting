classdef WfMatrix < mysort.wf.WfManager
    % This class manages detected events and cut multichannel waveform 
    % data. It does not store the raw data, but if a datasource is provided
    % that allows to access the raw data, recutting the waveforms is
    % possible.
    %
    % The class is conceptualized into
    %    single channel  DETECTION EVENTS
    %    multi  channel  WAVEFORMS
    %
    % Detection events of different channels can be grouped to form multi-
    % channel waveforms. Waveforms can be interpolated, subsample aligned,
    % subsample peak detected.

    properties         
        X
        cutLeft
        cutLength
        nC
    end
    
    methods
        %------------------------------------------------------------------
        function self = WfMatrix(X, eventTimes, eventChans, eventIDs, cutLeft, cutLength, ME)
            if nargin < 7
                ME = mysort.ds.MultiElectrode(unique(eventChans)', unique(eventChans));
            end
            self = self@mysort.wf.WfManager(eventTimes, eventChans, eventIDs, ME);
            self.X = X;
            self.cutLeft = cutLeft;
            self.cutLength = cutLength;
            self.nC = size(X,2)/cutLength;
            assert(size(X,2)/cutLength == round(size(X,2)/cutLength), 'cutLength does not match data matrix!');
        end
        %------------------------------------------------------------------
        function wfs = getWaveform4Idx(self, eventIdx, cutLeft, cutLength)
            wfs = self.X(eventIdx, :);
        end
        %------------------------------------------------------------------
        function wfManager = getWfManger4SubIdx(self, idx)
            wfManager = mysort.wf.WfMatrix(self.X(idx,:), self.eventTimes(idx),...
                self.eventChans(idx), self.eventIDs(idx), self.cutLeft, self.cutLength);
        end
        %------------------------------------------------------------------
        function nC = getNChannels(self)
            nC = self.nC;
        end
    end
end