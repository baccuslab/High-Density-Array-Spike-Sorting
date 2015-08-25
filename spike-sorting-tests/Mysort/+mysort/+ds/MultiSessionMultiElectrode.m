classdef MultiSessionMultiElectrode < mysort.ds.MultiElectrode
    properties
        % mandatory
        SessionIdxList
        MultiElectrodeList
        SingleSessionElectrodeIndices
        % internal
%         electrodeMESources
        electrodeCardinality
    end
    
    
    methods
        %------------------------------------------------------------------        
        function self = MultiSessionMultiElectrode(MElist, SessionIdxList)
            self = self@mysort.ds.MultiElectrode([], []);
            if nargin == 1 && isstruct(MElist)
                self.fromStruct(MElist);
                return
            end
            assert(isa(MElist, 'mysort.ds.MultiElectrode'), 'MElist must be of type MultiElectrode!');
            assert(length(MElist) == length(SessionIdxList), 'Every multielectrode must have a sessionidx!');
            
%             self.electrodeMESources = java.util.Hashtable();
            self.electrodeCardinality = [];
            for i=1:length(MElist)
                self.merge(MElist(i), SessionIdxList(i));
            end
        end
        %------------------------------------------------------------------        
        function fromStruct(self, S)
            self.SessionIdxList = S.SessionIdxList;
            self.SingleSessionElectrodeIndices = S.SingleSessionElectrodeIndices;
            self.electrodeCardinality = S.electrodeCardinality;
            self.electrodePositions = S.electrodePositions;
            self.electrodeNumbers = S.electrodeNumbers;
            self.electrodeLabels = S.electrodeLabels;
            self.parentElectrodeIndex = S.parenetElectrodeIndex;            
            for i=1:length(self.MultiElectrodeList)
                self.MultiElectrodeList(i) = mysort.ds.MultiElectrode(S.MultiElectrodeList(i));
            end
        end
        %------------------------------------------------------------------        
        function S = toStruct(self)
            S.SessionIdxList = self.SessionIdxList;
            S.SingleSessionElectrodeIndices = self.SingleSessionElectrodeIndices;
            S.electrodeCardinality = self.electrodeCardinality;
            S.electrodePositions = self.electrodePositions;
            S.electrodeNumbers = self.electrodeNumbers;
            S.electrodeLabels = self.electrodeLabels;
            S.parentElectrodeIndex = self.parentElectrodeIndex;
            for i=1:length(self.MultiElectrodeList)
                S.MultiElectrodeList(i) = self.MultiElectrodeList(i).toStruct();
            end
            S.readme = 'This struct was created by ds.MultiSessionMultiElectrode. Dont edit if you dont know what you are doing!';
        end 
        %------------------------------------------------------------------
        function save2File(self, fname, h5path)
            S = self.toStruct();
            mysort.h5.recursiveSave(fname, S, h5path);
        end
        %------------------------------------------------------------------
        function self = merge(self, ME, sessionIdx)
            assert(isa(ME, 'mysort.ds.MultiElectrode'), 'ME must be a multi electrode!');
            assert(isempty(self.SessionIdxList) || isempty(intersect(self.SessionIdxList,sessionIdx)), 'This session is already contained!');
            if ~isempty(self.electrodePositions)
                assert(~isempty(self.electrodeNumbers), 'You can only merge multielectrodes that have electrode numbers!');
            end            
            if ~isempty(ME.electrodePositions)
                assert(~isempty(ME.electrodeNumbers), 'You can only merge multielectrodes that have electrode numbers!');
            end

            bMergeLabels = true;
            if (~isempty(self.electrodePositions) && isempty(self.electrodeLabels)) ||...
               (~isempty(ME.electrodePositions)   && isempty(ME.electrodeLabels))
                bMergeLabels = false;
            end
            ElectrodeIndices = zeros(1, length(ME.electrodeNumbers));
            for i = 1:length(ME.electrodeNumbers)
                existingIdx = find(self.electrodeNumbers == ME.electrodeNumbers(i), 1);
                if isempty(existingIdx)
                    % New channel
                    ElectrodeIndices(i) = length(self.electrodeCardinality)+1;
                    self.electrodeCardinality = [self.electrodeCardinality 1];
                    self.electrodeNumbers = [self.electrodeNumbers ME.electrodeNumbers(i)];
                    self.electrodePositions  = [self.electrodePositions; ME.electrodePositions(i,:)];
                    if bMergeLabels
                        self.electrodeLabels = [self.electrodeLabels ME.electrodeLabel{i}];
                    end
                else
                    ElectrodeIndices(i) = existingIdx;
                    % Channel already existed
                    assert(all(self.electrodePositions(existingIdx,:) == ME.electrodePositions(i,:)), 'Tried to merge two electrodes with same Number but different positions!!!');
                    if bMergeLabels
                        assert(strcmp(self.electrodeLabels{existingIdx}, ME.electrodeLabel{i}), 'Tried to merge two electrodes with same Number but different Labels !!');
                    end
                    self.electrodeCardinality(existingIdx) = self.electrodeCardinality(existingIdx)+1;
                end
                % Add the current session to the electrodes session list
%                 key = num2str(ME.electrodeNumbers(i));
%                 list = sessionIdx;
%                 if self.electrodeMESources.containsKey(key)
%                     list = [self.electrodeMESources.get(key); list];
%                 end
%                 self.electrodeMESources.put(key, list);                
            end
            self.SessionIdxList     = [self.SessionIdxList sessionIdx];
            self.MultiElectrodeList = [self.MultiElectrodeList ME];
            self.SingleSessionElectrodeIndices = [self.SingleSessionElectrodeIndices {ElectrodeIndices}];
        end
        %------------------------------------------------------------------
        function c = getElectrodeCardinality(self, elnumber)
            if nargin == 1 || isempty(elnumber)
                c = self.electrodeCardinality;
            else
                elidx = self.getElIdx4ElNumber(elnumber);
                c = self.electrodeCardinality(elidx);
            end
        end
        %------------------------------------------------------------------
        function n = getNSessions(self)
            n = length(self.SessionIdxList);
        end
        %------------------------------------------------------------------
        function idx = getSessionsIdx(self)
            idx = self.SessionIdxList;
        end
        %------------------------------------------------------------------
        function nCs = getNChannelsPerSessions(self)
            n = self.getNSessions();
            nCs = zeros(1, n);
            for i=1:n
                nCs(i) = self.MultiElectrodeList(i).getNElectrodes();
            end
        end        
        %------------------------------------------------------------------
        function ME = getSubElectrode(self, elNumbers)
            nS = self.getNSessions();
            MElist = mysort.ds.MultiElectrode.empty();
            MEsessionidx = [];
            for i=1:nS
                MElist(end+1) = self.MultiElectrodeList(i).getSubElectrode(elNumbers);
                MEsessionidx(end+1) = self.SessionIdxList(i);
            end
            ME = mysort.ds.MultiSessionMultiElectrode(MElist, MEsessionidx);
        end  
        %------------------------------------------------------------------
        function ME = getSubSessionIdxMultiElectrode(self, global_sessionidx)
            local_sessionidx = self.globalSessionIdx2LocalSessionIdx(global_sessionidx);
            for i=1:length(local_sessionidx)
                MElist(i) = self.MultiElectrodeList(local_sessionidx(i));
            end
            ME = mysort.ds.MultiSessionMultiElectrode(MElist, global_sessionidx);
%             ME.parentElectrodeIndex = elidx;
        end   
        %------------------------------------------------------------------
        function lSIdx = globalSessionIdx2LocalSessionIdx(self, gSIdx)
            [a lSIdx c] = intersect(self.SessionIdxList, gSIdx);
        end
        %------------------------------------------------------------------
        function sIdx = getSessionIdx4ElectrodeNumber(self, elnumber)
            % returns the sessions in which all of the elnumbers are
            % present
            sIdx = find(all(cell2mat(arrayfun(@(y) {cellfun(@(x) any(ismember(self.electrodeNumbers(x), y)), self.SingleSessionElectrodeIndices)}, elnumber)')));
        end
%         %------------------------------------------------------------------
%         function lSIdx = getSessionIdx4ElectrodeIdx(self, elidx)
%             cellfun(@(x) any(ismember(ame.electrodeNumbers(x), elnumber)), ame.SingleSessionElectrodeIndices)
%         end
        
        
        %------------------ ELECTRODE POSITIONS ---------------------------
        %------------------------------------------------------------------
%         function setElectrodePositions(self, electrodePositions, electrodeNumbers, electrodeLabels)
%             error('not implemented');
%         end
%         %------------------ DATASOURCE ------------------------------------
%         function setDataSource(self, ds)
%             error('not implemented');
%         end    
%         %------------------ ELECTRODE LABELS ------------------------------
%         function setElectrodeLabels(self, el)
%             error('not implemented');
%         end
%         %------------------ ELECTRODE NUMBERS -----------------------------
%         function setElectrodeNumbers(self, en)
%             error('not implemented');
%         end        
    end
end