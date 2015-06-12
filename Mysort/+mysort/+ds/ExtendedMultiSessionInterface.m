classdef ExtendedMultiSessionInterface < mysort.ds.MultiSessionInterface & ...
                                         mysort.ds.ExtendedDataSourceInterface
            
    properties
        bufferedDetectedSpikes
    end
    
    methods(Abstract)

    end    
   
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = ExtendedMultiSessionInterface(name, s_per_sec, sessionList)
            assert(isa(sessionList, 'mysort.ds.DataSourceInterface'), 'sessionList must contain objects derived from type mysort.ds.DataSourceInterface!');
            if isempty(name)
                name = 'ExtendedMutliSessionInterface';
            end
            self = self@mysort.ds.MultiSessionInterface(name, s_per_sec, sessionList);
            self = self@mysort.ds.ExtendedDataSourceInterface(name, s_per_sec, sessionList(1).getMultiElectrode());
        end
        %------------------------------------------------------------------
        function setActiveSession(self, n)
            if isempty(self.activeSessionIdx) || isempty(n) || n~=self.activeSessionIdx
                self.bufferedDetectedSpikes = [];
                setActiveSession@mysort.ds.MultiSessionInterface(self, n);                
            end
        end
        %------------------------------------------------------------------
        function setSelectedSessions(self, sidx)
            setSelectedSessions@mysort.ds.MultiSessionInterface(self, sidx);
            self.bufferedDetectedSpikes = [];
        end   
       %------------------------------------------------------------------
        function [smad] = noiseStd(self, varargin)
            sessionIndex = self.activeSessionIdx;
            smad = self.sessionList(sessionIndex).noiseStd(varargin{:});
        end
        %------------------------------------------------------------------
        function [times pks] = detectSpikes(self, varargin)
            sessionIndex = self.activeSessionIdx;
            [times pks] = self.sessionList(sessionIndex).detectSpikes(varargin{:});   
        end
        %------------------------------------------------------------------
        function D = getSelectedSessionsDetectedSpikes(self)
            self.preLoadSelectedSessionsDetectedSpikes();
            D = self.bufferedDetectedSpikes;
        end
        %------------------------------------------------------------------
        function preLoadSelectedSessionsDetectedSpikes(self)
            if ~isempty(self.bufferedDetectedSpikes)
                return
            end
            sl = self.getSelectedSessions();
            D_ = [];
            for i=1:length(sl)
                [times pks] = sl(i).detectSpikes();
                ME = sl(i).getMultiElectrode;
                allTimes  = double(cell2mat(times));
                allTimes  = allTimes(:);
                allAmps   = double(cell2mat(pks));
                allAmps   = allAmps(:);
                allChans  = zeros(size(allTimes, 1), 1);
                nCum = 0;
                for k=1:length(times)
                    allChans(nCum+1:nCum+length(times{k}),1) = k;
                    nCum = nCum + length(times{k});
                end
                allChans = ME.electrodeNumbers(allChans);
                allTimes = allTimes(:)';
                allAmps = allAmps(:)';
                allChans = allChans(:)';
                idx = allTimes > 150 & allTimes < size(sl(i),1)-150;
                D_(i).allTimes = allTimes(idx);
                D_(i).allAmps = allAmps(idx);
                D_(i).allChans = allChans(idx); 
                D_(i).allSessions = ones(1, length(D_(i).allChans))*self.selectedSessions(i);
            end   
            D = struct;
            D.allTimes = [D_(:).allTimes];
            D.allAmps = [D_(:).allAmps];
            D.allChans = [D_(:).allChans];  
            D.allSessions = [D_(:).allSessions];  
            self.bufferedDetectedSpikes = D;
        end        
    end
end