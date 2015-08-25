classdef MultiSessionBufferedWfManager < mysort.wf.WfManagerWithDataSource
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
        eventSessions  % stores the sessionidx into the wfDataSource
                       % in which the spike was detected
        wfSessionIdx
        wfSingleSessionWfBuffers       % wfs per Session
        wfbufferCutLeft
        wfbufferCutLength
    end
    
    methods
        %------------------------------------------------------------------
        function self = MultiSessionBufferedWfManager(wfDataSource, eventTimes, eventChans, eventSessions, cutleft, cutlength, eventIDs, ME)
            if nargin < 8
                ME = [];
            end
            if nargin < 7 || isempty(eventIDs)
                eventIDs = 1:length(eventTimes);
            end        
            uniqSess = unique(eventSessions);
            self = self@mysort.wf.WfManagerWithDataSource(wfDataSource, eventTimes, eventChans, eventIDs, ME);
            self.wfSessionIdx = uniqSess;
            self.eventSessions = eventSessions;
            self.wfbufferCutLeft = cutleft;
            self.wfbufferCutLength = cutlength;
            self.wfSingleSessionWfBuffers = mysort.wf.BufferedWfManager.empty();
%             assert(~any(eventChans==0), 'Event channels must be greater zero!');
            for i=1:length(self.wfSessionIdx)
                idx = find(eventSessions == self.wfSessionIdx(i));
                self.wfSingleSessionWfBuffers(i) = mysort.wf.BufferedWfManager(...
                    wfDataSource.sessionList(self.wfSessionIdx(i)), eventTimes(idx), eventChans(idx), eventIDs(idx), cutleft, cutlength);
            end
            assert(~isempty(self.MultiElectrode), 'A MultiElectrode must be set!');
            assert(isa(self.MultiElectrode, 'mysort.ds.MultiSessionMultiElectrode'), 'The MultiElectrode must be a multiSession Electrode!');
        end
        
        %------------------------------------------------------------------
        function setCutLeftAndLength(self, cutleft, cutlength)
            for i=1:length(self.wfSingleSessionWfBuffers)
                wfb = self.wfSingleSessionWfBuffers(i);
                wfb.setCutLeftAndLength(cutleft, cutlength);
            end
        end
        %------------------------------------------------------------------
        function wfs = getWaveform4Idx(self, idx, varargin)
            if islogical(idx)
                idx = find(idx);
            end
            nS = length(self.wfSessionIdx);
            eS = self.eventSessions(idx);
            wfs = cell(nS,1);
            for i=1:nS
                sIdx = self.wfSessionIdx(i);
                B = self.wfSingleSessionWfBuffers(i);
                wfs{i} = B.getWaveform4ID(idx(eS==sIdx), varargin{:});
            end
        end
        %------------------------------------------------------------------
        function [T, ME, wfsAllSessions, nWfsPerElectrode] = getTemplate4Idx(self, idx, varargin)
            
            sessionIdx = self.wfSessionIdx;
            ME = self.MultiElectrode.getSubSessionIdxMultiElectrode(sessionIdx);
            if nargin == 5 && ~isempty(varargin{3})
                ME = ME.getSubElectrode(varargin{3});
                wfsAllSessions = self.getWaveform4Idx(idx, varargin{1:2}, ME.electrodeNumbers);
            else
                wfsAllSessions = self.getWaveform4Idx(idx, varargin{:});
            end
            
            
            nElectrodes = ME.getNElectrodes();
            elNumbers = ME.electrodeNumbers;
            nS = ME.getNSessions();
            
            T = zeros(nElectrodes, self.wfbufferCutLength);
            
            nWfs = zeros(1, nS);
            nWfsPerElectrode = zeros(1, nElectrodes);
            
            % The single session waveforms are in the order of the
            % electrodes of the multielectrode of the respective session.
            % Now we want to have a matrix with the templates with ALL
            % electrodes of ALL sessions. Thus we need a new multielectrode
            % with all those electrodes (we have that already, that is ME)
            % and fill the template matrix (T) in the appropriate order
            for i=1:length(sessionIdx)
                sIdx = sessionIdx(i);
                wfs = wfsAllSessions{i};
                nWfs(i) = size(wfs,1);
                if isempty(wfs)
                    continue
                end
                [session_elidx my_elidx] = ME.MultiElectrodeList(i).getElIdx4ElNumber(elNumbers); 
                
                if ~isempty(my_elidx)
                    nWfsPerElectrode(1, my_elidx) = nWfsPerElectrode(1, my_elidx) + nWfs(i);
                    T(my_elidx, :) = T(my_elidx, :) + mysort.wf.v2m(median(wfs, 1), length(my_elidx))*nWfs(i);                
                end
            end
            for i=1:nElectrodes
                T(i,:) = T(i,:)/nWfsPerElectrode(1,i);
            end
%             assert(all(~isnan(T(:))), 'nan in template calculation!');
        end
        %------------------------------------------------------------------
        function wfManager = getWfManger4SubIdx(self, idx)
            wfManager = mysort.wf.MultiSessionChildWfManager(self, idx);
%             wfManager = mysort.wf.MultiSessionBufferedWfManager(self.wfDataSource,...
%                 self.eventTimes(idx), self.eventChans(idx), self.eventSessions(idx), ...
%                 self.wfbufferCutLeft, self.wfbufferCutLength, self.eventIDs(idx));
%             nS = self.wfDataSource.getNSessions();
%             eS = self.eventSessions(idx);
%             for i=1:nS
%                 wfManager.wfSingleSessionWfBuffers(i) = self.wfSingleSessionWfBuffers(i).getWfManger4SubIdx(idx(eS==i));
%             end            
        end          
    end
end