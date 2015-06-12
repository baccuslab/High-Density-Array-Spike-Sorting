classdef MultiElectrode < handle
    properties
        % internal 
        parentElectrodeIndex
        
        % mandatory
        electrodePositions
        electrodeNumbers
        
        % optional        
        dataSource
        electrodeLabels        
    end
    
    
    methods
        %------------------------------------------------------------------        
        function self = MultiElectrode(electrodePositions, electrodeNumbers, varargin)
            if nargin == 1 && isstruct(electrodePositions)
                self.fromStruct(electrodePositions);
                return
            end
            P.dataSource = [];
            P.electrodeLabels = [];
            P.parentElectrodeIndex = [];
            P = mysort.util.parseInputs(P, varargin, 'error');
            self.setElectrodePositions(electrodePositions, electrodeNumbers);
            
            if ~isempty(P.dataSource)
                self.setDataSource(P.dataSource);
            end
            self.electrodeLabels = P.electrodeLabels;
            self.parentElectrodeIndex = P.parentElectrodeIndex;
        end
        %------------------------------------------------------------------        
        function fromStruct(self, S)
            self.setElectrodePositions(S.electrodePositions, S.electrodeNumbers);
            self.electrodeLabels = S.electrodeLabels;
            self.parentElectrodeIndex = S.parentElectrodeIndex;
        end
        %------------------------------------------------------------------        
        function S = toStruct(self)
            S.electrodePositions = self.electrodePositions;
            S.electrodeNumbers = self.electrodeNumbers;
            S.electrodeLabels = self.electrodeLabels;
            S.parentElectrodeIndex = self.parentElectrodeIndex;
        end 
        %------------------------------------------------------------------
        function save2File(self, fname, h5path)
            S = self.toStruct();
            mysort.h5.recursiveSave(fname, S, h5path);
        end
        %------------------------------------------------------------------
        function cp = copy(self)
            S = self.toStruct();
            cp = mysort.ds.MultiElectrode(S);
        end
        
        %------------------------------------------------------------------
        function [merged mergedElIdx] = merge(self, ME)
            % Merges this m-electrode with the m-electrode in ME.
            % mergedElIdx is the array that contains for each electrode in
            % the resulting merged electrode the index into the old
            % electrode set where it came from or the new electrode idx if 
            % this electrode was only in ME
            merged = self.copy();
            assert(isa(ME, 'mysort.ds.MultiElectrode'), 'ME must be a multi electrode!');
            if ~isempty(self.electrodePositions)
                assert(~isempty(self.electrodeNumbers), 'You can only merge multielectrodes that have electrode numbers!');
            end
            if ~isempty(ME.electrodePositions)
                assert(~isempty(ME.electrodeNumbers), 'You can only merge multielectrodes that have electrode numbers!');
            end
            merged.dataSource = [];
            bMergeLabels = true;
            if (~isempty(self.electrodePositions) && isempty(self.electrodeLabels)) ||...
               (~isempty(ME.electrodePositions)   && isempty(ME.electrodeLabels))
                bMergeLabels = false;
            end
            mergedElIdx = zeros(ME.getNElectrodes(),1);
            for i=1:length(ME.electrodeNumbers)
                myElIdx = find(merged.electrodeNumbers == ME.electrodeNumbers(i));
                if isempty(myElIdx)
                    mergedElIdx(i) = length(merged.electrodeNumbers)+1;
                    
                    merged.electrodeNumbers = [merged.electrodeNumbers(:); ME.electrodeNumbers(i)];
                    merged.electrodePositions = [merged.electrodePositions;
                                                 ME.electrodePositions(i,:)];
                    if bMergeLabels
                        merged.electrodeLabels = [merged.electrodeLabels ME.electrodeLabels{i}];
                    end
                else
                    assert(length(myElIdx) == 1, 'More than one electrode found!');
                    mergedElIdx(i) = myElIdx;
                end
            end
        end
        %------------------------------------------------------------------
        function ME = getSubElectrode(self, elnumbers)
            ME = self.getSubElectrode4ElIdx(self.getElIdx4ElNumber(elnumbers));
        end 
        %------------------------------------------------------------------
        function ME = getSubElectrode4ElIdx(self, elidx)
            epos = self.electrodePositions(elidx,:);
            elab = [];
            if ~isempty(self.electrodeLabels)
                elab = self.electrodeLabels(elidx);
            end
            enum = self.electrodeNumbers(elidx);
            epar = elidx;
            if ~isempty(self.parentElectrodeIndex)
                epar = self.parentElectrodeIndex(elidx);
            end
            ME = mysort.ds.MultiElectrode(epos, enum, ...
                'dataSource', [], ...
                'electrodeLabels', elab, ...
                'parentElectrodeIndex', epar);
        end        
        %------------------------------------------------------------------
        function [elidx calleridx found_elNumbers] = getElIdx4ElNumber(self, elNumbers)
            % This function returns for a given set of electrode numbers
            % the indices of those numbers into the own electrode numbers (elidx),
            % but also the indices into elNumbers of the electrodes that are 
            % actually in this electrode
            elidx = -1*ones(1, length(elNumbers));
            for i=1:length(elNumbers)
                el_idx_in_mylist = find(self.electrodeNumbers == elNumbers(i), 1);
                if ~isempty(el_idx_in_mylist)
                    elidx(i) = el_idx_in_mylist;
                end
            end
            calleridx = []; count = 1;
            for i=1:length(self.electrodeNumbers)
                el_idx_in_callerlist = find(self.electrodeNumbers(i) == elNumbers, 1);
                if ~isempty(el_idx_in_callerlist)
                    calleridx(count) = el_idx_in_callerlist; 
                    count = count+1;
                end
            end

            found_elNumbers = elNumbers(elidx>0);
            elidx = elidx(elidx>0);
            %[elnum elidx caller_elidx] = intersect(self.electrodeNumbers, elnumbers);
        end
       
        %------------------ ELECTRODE POSITIONS ---------------------------
        %------------------------------------------------------------------
        function setElectrodePositions(self, electrodePositions, electrodeNumbers, electrodeLabels)
            if isempty(electrodePositions)
                return
            end
            assert(isnumeric(electrodePositions), 'electrodePositions must be numeric!');
            assert(any([1 2 3] == size(electrodePositions,2)), ...
                'Electrode positions must be one, two or three dimensional!');
            self.electrodePositions = electrodePositions;
            self.setElectrodeNumbers(electrodeNumbers);
            if nargin > 3
                self.setElectrodeLabels(electrodeLabels);
            else
                self.setElectrodeLabels([]);
            end
        end
        %------------------------------------------------------------------
        function epos = getElectrodePositions(self, elnumber)
            if nargin == 1
                epos = self.getElectrodePositions4ElIdx();
                return
            end            
            epos = self.getElectrodePositions4ElIdx(self.getElIdx4ElNumber(elnumber));
        end
        %------------------------------------------------------------------
        function epos = getElectrodePositions4ElIdx(self, elidx)
            if nargin == 1
                elidx = 1:self.getNElectrodes();
            end
            epos = self.electrodePositions(elidx,:);
        end        
        %------------------ DATASOURCE ------------------------------------
        %------------------------------------------------------------------
        function setDataSource(self, ds)
            if isempty(ds)
                self.dataSource = [];
                return
            end
%             assert(self.getNElectrodes() == ds.getNChannels(),'Every data channel must have an electrode position!');            
            self.dataSource = ds;
        end
        %------------------------------------------------------------------
        function ds = getDataSource(self)
            ds = self.dataSource;
        end    
        %------------------------------------------------------------------
        function b = hasDataSource(self)
            b = ~isempty(self.dataSource);
        end        
        
        %------------------ ELECTRODE LABELS ------------------------------
        %------------------------------------------------------------------
        function el = getElectrodeLabels(self)
            el = self.electrodeLabels;
        end
        %------------------------------------------------------------------
        function setElectrodeLabels(self, el)
            if isempty(el)
                self.electrodeLabels = el;
                return
            end
            assert(self.getNElectrodes() == length(el),'Every data channel must have an label!');            
            self.electrodeLabels = el;
        end
        %------------------------------------------------------------------
        function b = hasElectrodeLabels(self)
            b = ~isempty(self.electrodeLabels);
        end        
        %------------------ ELECTRODE NUMBERS -----------------------------
        %------------------------------------------------------------------
        function en = getElectrodeNumbers(self)
            en = self.electrodeNumbers;
        end
        %------------------------------------------------------------------
        function setElectrodeNumbers(self, en)
            if isempty(en)
                self.electrodeNumbers = en;
                return
            end
            assert(self.getNElectrodes() == length(en),'Every data channel must have a number!');            
            self.electrodeNumbers = en;
        end        
        %------------------------------------------------------------------
        function b = hasElectrodeNumbers(self)
            b = ~isempty(self.electrodeNumbers);
        end       


        %------------------ UTILITY FUNCTIONS -----------------------------
        %------------------------------------------------------------------
        function b = eq(A, B)
            b = true;
            [nElA nDimA] = size(A.electrodePositions);
            [nElB nDimB] = size(B.electrodePositions);
            if nElA ~= nElB || nDimA ~=nDimB
                b = false;
                return
            end
            if any(A.electrodePositions(:) ~= B.electrodePositions(:))
                b = false;
                return
            end
            if length(A.electrodeNumbers) ~= length(B.electrodeNumbers)
                b = false;
                return
            end
            if ~isempty(A.electrodeNumbers) && any(A.electrodeNumbers ~= B.electrodeNumbers)
                b = false;
                return
            end
        end
        
        function n = getNElectrodes(self)
            n = size(self.electrodePositions,1);
        end
        %------------------------------------------------------------------        
        function n = getNDimensions(self)
            n = size(self.electrodePositions,2);
        end

        %------------------------------------------------------------------        
        function d = getDistance(self, e1, e2)
            d = norm(self.electrodePositions(e1,:) ...
                    -self.electrodePositions(e2,:));
        end
        %------------------------------------------------------------------        
        function d = getAllElDistancesFromPoint(self, point)
            if isempty(point)
                d = [];
                return
            end
            d = sqrt(sum((repmat(point,self.getNElectrodes(),1) - self.electrodePositions).^2,2)); 
        end
        %------------------------------------------------------------------        
        function N = getNeighbors(self, e1, maxDist)
            assert(~isempty(e1), 'e1 must not be empty+');
            if nargin == 2
                maxDist = 20;
            end
            myPos = self.electrodePositions(e1,:);
            dists = self.getAllElDistancesFromPoint(myPos);
            N = find(dists<=maxDist);
        end 
        %------------------------------------------------------------------        
        function [groupsidx nGroupsPerElectrode] = getLocalElectrodeGroups(self)
            [groupsidx nGroupsPerElectrode] = mysort.mea.constructLocalElectrodeGroups(self.electrodePositions(:,1), self.electrodePositions(:,2));       
        end
        %------------------------------------------------------------------        
        function ah = plotElectrodeGroups(self, groupsidx)
            if nargin == 1
                [groupsidx nGroupsPerElectrode] = self.getLocalElectrodeGroups();
            end
            figure;
            ah = axes();
            plot(self.electrodePositions(:,1), self.electrodePositions(:,2), 'ok');
            hold on
            for i=1:length(groupsidx)
                x = self.electrodePositions(groupsidx{i},1);
                y = self.electrodePositions(groupsidx{i},2);
                plot(ah1, x+1+2*rand, y+2*rand, 'x', 'color', mysort.plot.vectorColor(i), 'markersize', 14, 'linewidth', 2);
                hold on
            end    
        end
    end
end