classdef MultiSessionInterface < mysort.ds.DataSourceInterface
    properties
        sessionList
        activeSessionIdx
        selectedSessions
        allSessionsMergedMultiElectrode
        allSessionLengths
        size_buffer
    end
    
    properties (Access=private)
        bIsConcatenated
    end
    
   
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = MultiSessionInterface(name, s_per_sec, sessionList)
            assert(isa(sessionList, 'mysort.ds.DataSourceInterface'), 'sessionList must contain objects derived from type mysort.ds.DataSourceInterface!');
            if isempty(name)
                name = 'MutliSessionInterface';
            end
            self = self@mysort.ds.DataSourceInterface(name, s_per_sec, sessionList(1).getMultiElectrode());
            self.sessionList = sessionList; 
            self.allSessionLengths = [];
            self.size_buffer = [];
            self.bIsConcatenated = false;
%             for i=1:length(sessionList)
%                 
%             end
            self.setActiveSession(1);
            self.concatenateSessions();
            % call this function once to make sure it is initialized after
            % the constructer is left. Had problems with Derk accessing the
            % variable before it was set.
            self.getAllSessionsLength();
        end
        

        
        %% TOP LEVEL DATA ACCESS FUNCTIONS     
        %------------------------------------------------------------------
        function out = end(self,k,n)
            if k<3
                out = self.size(k);
            elseif k==3
                out = self.getNSessions();
            else
                error('not implemented');
            end
        end 
        %------------------------------------------------------------------
        function varargout = size(self,varargin)
            if nargin == 2 && varargin{1} == 3
                varargout{1} = length(self.sessionList);
                return
            end

            if isempty(self.activeSessionIdx)
                dims = [0 0];
            elseif ~isempty(self.size_buffer)
                dims = self.size_buffer;
            else
                if ~self.bIsConcatenated
                    [rows cols] = self.sessionList(self.getActiveSessionIdx()).size();
                    dims = [rows cols self.getNSessions()];
                else
                    rows = sum(self.getAllSessionsLength());
                    cols = size(self.sessionList(1),2);
                    dims = [rows cols];
                end
            end
            varargout=matlabfilecentral.parseSize.parseSize(dims,nargout,varargin{:}); 
        end
        %------------------------------------------------------------------
        function [rows cols] = size_(self, dim, session_idx)
%             error('Is this function still used? 02.10.2013 Felix');
            if nargin == 1
                dim = [];
                session_idx = self.activeSessionIdx;
            elseif nargin == 2
                session_idx = self.activeSessionIdx;
            end
            
            [rows cols] = self.sessionList(session_idx).size();
            
            if nargout == 0 && isempty(dim)
                rows = [rows cols];
            elseif nargout <= 1 && dim == 2
                rows = cols;
            end              
        end
        %------------------------------------------------------------------
        function L = getNSamples_(self, dim, session_idx)
            if nargin == 1
                dim = [];
                session_idx = self.activeSessionIdx;
            elseif nargin == 2
                session_idx = self.activeSessionIdx;
            end
            if ~self.bIsConcatenated
                [rows cols] = self.sessionList(session_idx).size();
            else
                rows = sum(structfun(@(x) x.size(1), self.sessionList));
            end
            
            L = rows;         
        end
        %------------------------------------------------------------------
        function s = getSessionList(self)
            s = self.sessionList;
        end
        %------------------------------------------------------------------
        function n = getNSessions(self)
%             if ~self.bIsConcatenated
                n = length(self.sessionList);
%             else
%                 n = 1;
%             end
        end
        %------------------------------------------------------------------
        function n = getNSelectedSessions(self)
            n = length(self.selectedSessions);
        end
        %------------------------------------------------------------------
        function idx = getSessionsIdx4Filenames(self, filenames)
            myFnames = self.getSessionFilenames();
            idx = zeros(1, length(filenames));
            for i=1:length(filenames)
                funh = @(x) ~isempty(strfind(filenames{i},x));
                idx(i) = find(cellfun(funh, myFnames));
            end
        end        
        %------------------------------------------------------------------
        function restrictToChannels(self, varargin)
            self.size_buffer = [];
            restrictToChannels@mysort.ds.DataSourceInterface(self, varargin{:});
            for i=1:length(self.sessionList)
                self.sessionList(i).restrictToChannels(varargin{:});
            end
        end   
        %------------------------------------------------------------------
        function concatenateSessions(self)
            self.size_buffer = [];
            for i=2:length(self.sessionList)
                assert(self.sessionList(1).MultiElectrode == self.sessionList(i).MultiElectrode, 'Sessions cannot be concatenated since they dont share the same MultiElectrode!');
            end
            self.bIsConcatenated = true;
        end
        
        %------------------------------------------------------------------
        function s = getActiveSession(self)
            s = [];
            if isempty(self.sessionList)
                return
            end
            s = self.sessionList(self.activeSessionIdx);
        end
        %------------------------------------------------------------------
        function n = getActiveSessionIdx(self)
            n = self.activeSessionIdx;
        end
        %------------------------------------------------------------------
        function setActiveSession(self, n)
            assert(length(n) == 1, 'only one session can be active!')
            assert(n<=self.getNSessions(), 'not an allowed session index!')
            assert(n>0, 'not an allowed session index!')
            self.activeSessionIdx = n;
            self.MultiElectrode = self.sessionList(n).getMultiElectrode();
            if ~any(self.selectedSessions == n)
                self.selectedSessions = sort([self.selectedSessions n]);
            end
        end
        %------------------------------------------------------------------
        function sidx = getSelectedSessionIdx(self)
            sidx = self.selectedSessions;
        end
        %------------------------------------------------------------------
        function sl = getSelectedSessions(self)
            sl = self.sessionList(self.selectedSessions);
        end
        %------------------------------------------------------------------
        function me = getSelectedSessionsMergedMultiElectrode(self)
            sidx = self.getSelectedSessionIdx();
            mes = self.getSelectedSessionsMultiElectrodes();
            me = mysort.ds.MultiSessionMultiElectrode(mes, sidx);
            me.dataSource = self;
        end
        %------------------------------------------------------------------
        function mes = getSelectedSessionsMultiElectrodes(self)
            i = self.getSelectedSessionIterator();
            mes = mysort.ds.MultiElectrode.empty();
            while i.hasNext()
                s = i.next();
                mes(i.idx) = s.getMultiElectrode();
            end
        end        
        %------------------------------------------------------------------
        function me = getMergedMultiElectrode4Sessions(self, sidx)
            mes = self.getMultiElectrodes4Sessions(sidx);
            me = mysort.ds.MultiSessionMultiElectrode(mes, sidx);
            me.dataSource = self;
        end           
        %------------------------------------------------------------------
        function mes = getMultiElectrodes4Sessions(self, sidx)
            i = mysort.util.Iterator(self.sessionList(sidx));
            mes = mysort.ds.MultiElectrode.empty();
            while i.hasNext()
                s = i.next();
                mes(i.idx) = s.getMultiElectrode();
            end
        end         
        %------------------------------------------------------------------
        function offsets = getSelectedSessionsTimeOffsets(self)
            L = self.getSelectedSessionLength();
            offsets = cumsum(L);
            offsets = [0 offsets(1:end-1)];            
        end
        %------------------------------------------------------------------
        function L = getSelectedSessionLength(self)
            L = zeros(1,self.getNSessions());
            i = self.getSelectedSessionIterator();
            while i.hasNext()
                s = i.next(); 
                L(i.idx) = size(s,1);
            end
        end
        %------------------------------------------------------------------
        function L = getAllSessionsLength(self)
            if ~isempty(self.allSessionLengths)
                L = self.allSessionLengths;
                return
            end
            L = zeros(1, length(self.sessionList));
            for i=1:length(self.sessionList)
                L(i) = self.sessionList(i).size(1);
            end
            self.allSessionLengths = L;
        end   
        %------------------------------------------------------------------
        function mes = getAllSessionsMultiElectrodes(self)
            mes = mysort.ds.MultiElectrode.empty();
            for i=1:length(self.sessionList)
                mes(i) = self.sessionList(i).getMultiElectrode();
            end
        end            
        %------------------------------------------------------------------
        function me = getAllSessionsMergedMultiElectrode(self)
            if ~isempty(self.allSessionsMergedMultiElectrode)
                me = self.allSessionsMergedMultiElectrode;
            else
                me = mysort.ds.MultiSessionMultiElectrode(self.sessionList(1).getMultiElectrode(), 1);
                for i=2:length(self.sessionList)
                    me.merge(self.sessionList(i).getMultiElectrode(), i);
                end
                me.dataSource = self;
                self.allSessionsMergedMultiElectrode = me;
            end
        end        
                    
        %------------------------------------------------------------------
        function setSelectedSessions(self, sidx)
            assert(~isempty(sidx), 'sidx must not be empty!')
            assert(all(sidx > 0), 'session index cant be zero or negative!')
            assert(all(sidx <= self.getNSessions()), 'Session index out of bounds!')
            self.selectedSessions = sidx;
            self.setActiveSession(sidx(1));
        end         
        
        %------------------------------------------------------------------
        function erg = getSessionVar(self, varname, sessionidx)
            if nargin == 2
                sessionidx = 1:self.getNSessions();
            end
            sessions = self.sessionList(sessionidx);
            if ischar(sessions(1).(varname))
                erg = cell(1, length(sessions));
                for i=1:length(sessions)
                    erg{i} = sessions(i).(varname);
                end
            else
                erg = zeros(1, length(sessions));
                for i=1:length(sessions)
                    erg(i) = sessions(i).(varname);
                end
            end
        end
        %------------------------------------------------------------------
        function L = getSessionLength(self, sidx)
            L = zeros(1,length(sidx));
            i = self.getSessionIterator(sidx);
            while i.hasNext()
                s = i.next(); 
                L(i.idx) = size(s,1);
            end
        end        
        %------------------------------------------------------------------
        function Session = getSessionObject(self, sessionIndex)
            if nargin == 1
                sessionIndex = 1;
            end
            Session = self.sessionList(sessionIndex);
        end        
        %------------------------------------------------------------------
        function iter = getSessionIterator(self, sessionIndex)
            if nargin == 1
                sessionIndex = 1:self.size(3);
            end
            assert(isnumeric(sessionIndex), 'SessionIndex must be numeric!');
            iter = mysort.util.Iterator(self.sessionList(sessionIndex));
        end    
        %------------------------------------------------------------------
        function iter = getSelectedSessionIterator(self)
            sl = self.getSelectedSessions();
            iter = mysort.util.Iterator(sl);
        end         
        
        %------------------------------------------------------------------
        function wf = getWaveform(self, t, cutLeft, cutLength, channels, sessionIdx)
            if nargin < 5 || isempty(channels)
                channels = 1:self.size_(2);
            end
            if ~self.bIsConcatenated
                if nargin < 6
                    sessionIdx = self.activeSessionIdx;
                else
                    assert(sessionIdx > 0 & sessionIdx <= length(self.sessionList), sprintf('sessionIdx out of bounds! %d', sessionIdx));
                end
                wf = self.sessionList(sessionIdx).getWaveform(t, cutLeft, cutLength, channels);
            else
                wf = self.getWaveformConcatenated_(t, cutLeft, cutLength, channels);
            end
        end

        %------------------------------------------------------------------
        function S = getSpikeSorting(self, sorting_idx, sessionidx)
            if isempty(sorting_idx)
                warning('requested empty sorting_idx');
                return
            end
            if nargin == 2 || isempty(sessionidx)
                sessionidx = self.activeSessionIdx;
            end
            assert(length(sessionidx) == 1, 'only one session can be requested for single session spike sortings');
            S = self.SpikeSortingContainers{sorting_idx}.getSingleSessionSorting(sessionidx, self.MultiElectrode);
            S.wfDataSource = self;
        end
        %------------------------------------------------------------------
        function S = getMultiSessoinSpikeSorting(self, sorting_idx)
            S = self.SpikeSortingContainers{sorting_idx};
        end

    end
    methods(Sealed)
        %------------------------------------------------------------------
        function X = getData(self, varargin)
            if ~self.bIsConcatenated || nargin == 4
                X = self.getData_(varargin{:});
            else
                X = self.getDataConcatenated_(varargin{:});
            end
        end
        %------------------------------------------------------------------
        function X = getData_(self, timeIndex, channelIndex, sessionIndex)
            if nargin == 4 && length(sessionIndex)>1
                % concatenate data from multiple sessions!
                % TODO! allocate memory first !!!
                % TODO! check that all subsession indices are valid !
                X = self.getData_(timeIndex, channelIndex, sessionIndex(1));
                for i=2:length(sessionIndex)
                    X = [X; self.getData_(timeIndex, channelIndex, sessionIndex(i))];
                end
                return
            end
              
            if nargin < 4
                sessionIndex = self.activeSessionIdx;
                if nargin < 3
                    channelIndex = 1:self.size_(2);
                    if nargin < 2
                        timeIndex = 1:self.size_(1);
                    end
                end
            end
            
            if strcmp(timeIndex, ':')
                timeIndex = 1:self.sessionList(sessionIndex).size(1);
            end                    
            if strcmp(channelIndex, ':')
                channelIndex = 1:self.sessionList(sessionIndex).size(2);
            end                    

            X = self.sessionList(sessionIndex).getData(timeIndex, channelIndex);
        end
        %------------------------------------------------------------------
        function X = getDataConcatenated_(self, timeIndex, channelIndex) 
            if nargin < 3 || strcmp(channelIndex, ':')
                channelIndex = 1:self.size(2);            
            end          
            SL = self.getAllSessionsLength();
            if strcmp(timeIndex, ':')
                timeIndex = 1:sum(SL);
            end                    
         
            sessionOffsets = cumsum([0 SL]);
            X = zeros(length(timeIndex), length(channelIndex));
            totalLength = 0;
            for i=1:length(sessionOffsets)-1
                sessionFirstIndex = sessionOffsets(i)+1;
                sessionLastIndex  = sessionOffsets(i+1);
                timeIndexInThisSession = timeIndex(timeIndex>=sessionFirstIndex & timeIndex <= sessionLastIndex) - sessionFirstIndex +1;
                if ~isempty(timeIndexInThisSession)
%                     assert(min(timeIndexInThisSession) > 0, 'Invalid Time Index! (<1)');
%                     assert(max(timeIndexInThisSession) <= size(self.sessionList(i),1), 'Invalid Time Index! (>end)');                    
                    tmp = self.getData_(timeIndexInThisSession, channelIndex, i);
                    X(totalLength+1:totalLength+length(timeIndexInThisSession),:) = tmp;
                    totalLength = totalLength+length(timeIndexInThisSession);
                end
            end
        end
        %------------------------------------------------------------------
        function wfs = getWaveformConcatenated_(self, t, cutLeft, cutLength, channelIndex) 
            if nargin < 5 || isempty(channelIndex)
                channelIndex = 1:self.size(2);            
            end          
            SL = self.getAllSessionsLength();       
            sessionOffsets = cumsum([0 SL]);
       
            wfs = zeros(length(t), cutLength*length(channelIndex));

            for i=1:length(sessionOffsets)-1
                sessionFirstIndex = sessionOffsets(i)+1;
                sessionLastIndex  = sessionOffsets(i+1);
                inThisSessionIdx = t>=sessionFirstIndex & t <= sessionLastIndex;
                tInThisSession = t(inThisSessionIdx) - sessionFirstIndex +1;
                if ~isempty(tInThisSession)         
                    wfs_tmp = self.sessionList(i).getWaveform(tInThisSession, cutLeft, cutLength, channelIndex);
                    wfs(inThisSessionIdx, :) = wfs_tmp;
                end
            end
        end        
    end
end