classdef BufferedWfManager < mysort.wf.WfManagerWithDataSource
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
        wfbuffer       % wfs that were already cut
        wfbufferIdx
        wfbufferN
        wfbufferCutLeft
        wfbufferCutLength
    end
    
    methods
        %------------------------------------------------------------------
        function self = BufferedWfManager(wfDataSource, eventTimes, eventChans, eventIDs, cutleft, cutlength)
            self = self@mysort.wf.WfManagerWithDataSource(wfDataSource, eventTimes, eventChans, eventIDs);
            self.setCutLeftAndLength(cutleft, cutlength);
        end
        
        %------------------------------------------------------------------
        function setCutLeftAndLength(self, cutleft, cutlength)
            % TODO: If new parameters are smaller than old, keep buffer
            self.wfbufferCutLeft = cutleft;
            self.wfbufferCutLength = cutlength;
            self.bufferReset();
        end
        %------------------------------------------------------------------
        function bufferReset(self)
            nC = self.wfDataSource.MultiElectrode.getNElectrodes();
            if isempty(self.eventTimes)
                self.wfbuffer = [];
                self.wfbufferIdx = [];
            else
                self.wfbuffer = zeros(500, nC*self.wfbufferCutLength);
                self.wfbufferIdx = sparse(length(self.eventTimes), 1, 0);
            end
            self.wfbufferN = 0;            
        end   
        %------------------------------------------------------------------
        function wfs = getWaveform4Idx(self, idx, cutleft, cutlength, electrodeNumber)
            if islogical(idx)
                idx = find(idx);
            end
            if isempty(idx)
                wfs = [];
                return
            end
%             assert(any(self.MultiElectrode.getElIdx4ElNumber(electrodeNumber)==-1), 'could not find requested electrode number!')
            bufferIdx = self.wfbufferIdx(idx,1);
            % Check which waveforms were not buffered yet and buffer
            notBufferedIdx = find(bufferIdx==0);
            % Check if buffer is big enough to hold all new waveforms
            if size(self.wfbuffer,1)-self.wfbufferN < length(notBufferedIdx)
                self.wfbuffer(end+max(500,length(notBufferedIdx)),1) = 0;
            end            
            if ~isempty(notBufferedIdx)
                fprintf('Buffering %d waveforms...', length(notBufferedIdx));
                tic
                wfIdx = idx(notBufferedIdx);
                bIdx  = self.wfbufferN + (1:length(notBufferedIdx));
                self.wfbufferIdx(wfIdx,1) = bIdx;
                self.wfbufferN = bIdx(end);
                wfs = getWaveform4Idx@mysort.wf.WfManagerWithDataSource(self,...
                        wfIdx, self.wfbufferCutLeft, self.wfbufferCutLength);
                %wfs = getWaveform4Idx@mysort.wf.WfManagerWithDataSource(self,wfIdx, self.wfbufferCutLeft, self.wfbufferCutLength);
                self.wfbuffer(bIdx,:) = wfs;
                tElapsed=toc;
                fprintf('%.4fsec - %.4fms per wf\n', tElapsed, 1000*tElapsed/length(notBufferedIdx));
            end
            % Retrieve waveforms from buffer
            if nargin < 5 || isempty(electrodeNumber)
                wfs = self.wfbuffer(self.wfbufferIdx(idx,1),:);
            else
                [channelidx] = self.MultiElectrode.getElIdx4ElNumber(electrodeNumber);
                chansubidx = mysort.wf.vSubChannelIdx(self.wfbufferCutLength, self.MultiElectrode.getNElectrodes(), channelidx);
                wfs = self.wfbuffer(self.wfbufferIdx(idx,1), chansubidx);
            end
        end   
        %------------------------------------------------------------------
        function wfManager = getWfManger4SubIdx(self, idx)
            wfManager = mysort.wf.BufferedWfManager(self.wfDataSource,...
                self.eventTimes(idx), self.eventChans(idx), self.eventIDs(idx),...
                self.wfbufferCutLeft, self.wfbufferCutLength);
            if ~isempty(idx)
                bIdx = self.wfbufferIdx(idx,1);
                alreadyBufferedIdx = bIdx==1;
                wfManager.wfbuffer = self.wfbuffer(bIdx(alreadyBufferedIdx),:);
                wfManager.wfbufferIdx = sparse(length(idx), 1, 0);
                wfManager.wfbufferN = sum(alreadyBufferedIdx);
                wfManager.wfbufferIdx(1:wfManager.wfbufferN) = find(alreadyBufferedIdx);
            end
        end        
    end
end