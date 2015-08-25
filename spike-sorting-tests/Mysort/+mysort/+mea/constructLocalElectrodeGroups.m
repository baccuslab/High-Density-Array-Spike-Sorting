function [groups nGroupsPerElectrode] = constructLocalElectrodeGroups(x, y, varargin)
    P.maxElPerGroup = 7;
    P.minElPerGroup = 1;
    P.addIfNearerThan = 20; % always add direct neighbors
%     P.maxOverlapPerGroup = 6;
    P.maxDistanceWithinGroup = 52;  %keep over 51.59 !
    P = mysort.util.parseInputs(P, varargin, 'error');

    groups = {};
    
    N = length(x);
    nGroupsPerElectrode = zeros(N,1);
    currentGroup = 1;
    while 1
        notInAnyGroup = nGroupsPerElectrode==0;
        if ~any(notInAnyGroup)
            break
        end
        notInAnyGroupIdx = find(notInAnyGroup);        
        % take the top left most electrode not in any group
        [s sidx] = sortrows([x(notInAnyGroupIdx) y(notInAnyGroupIdx)]);
        seed = notInAnyGroupIdx(sidx(1));
        % compute Distances
        D = computeDists([x(seed) y(seed)], [x y]);
        
        % first, add all electrodes, no matter what, that are closer than a
        % certain distance (direct neighbors)
        group = find(D<=P.addIfNearerThan);
        group = sortrows([group D(group)], 2);
        groups{currentGroup} = group(1:min(size(group,1),P.maxElPerGroup),1);
        
        % check if that is enough
        if size(groups{currentGroup},1) < P.maxElPerGroup
            % no, not enough
            % add nearest electrodes that are not in any group first
            D = max(computeDists([x(groups{currentGroup}) y(groups{currentGroup})], [x y]),[],2);
            group2 = find(D>P.addIfNearerThan & D<P.maxDistanceWithinGroup & notInAnyGroup);
            group2 = sortrows([group2 D(group2)], 2);
        
            if ~isempty(group2)
                groups{currentGroup} = [groups{currentGroup}; group2(1:min(end, P.maxElPerGroup-size(groups{currentGroup},1)),1)];
                D = max(computeDists([x(groups{currentGroup}) y(groups{currentGroup})], [x y]),[],2);
            end
            % check if that is enough
            if size(groups{currentGroup},1) < P.maxElPerGroup
                % now consider also electrodes that are already in
                % a group, this does not contain the seed !
                group3 = find(D>P.addIfNearerThan & D<P.maxDistanceWithinGroup & ~notInAnyGroup);
                group3 = sortrows([group3 D(group3)], 2);
                % check if this is enough
                assert(length(group) + length(group2) + length(group3) >= P.minElPerGroup,'cannot construct groups with these parameters!');

                groups{currentGroup} = [groups{currentGroup}; group3(1:min(end, P.maxElPerGroup-size(groups{currentGroup},1)),1)];
            end
        end
        nGroupsPerElectrode(groups{currentGroup}) = nGroupsPerElectrode(groups{currentGroup}) +1;            
        currentGroup = currentGroup + 1;           
    end
    
    % check if there are groups that are completely covered by other
    % groups. Those we can remove.
    removeCoveredGroups();
    mergeGroupsWithOnly1Difference();
    
    %----------------------------------------------------------------------
    function removeCoveredGroups()
        for i=length(groups):-1:1
            if all(nGroupsPerElectrode(groups{i})>1)
                nGroupsPerElectrode(groups{i}) = nGroupsPerElectrode(groups{i})-1;
                groups(i) = [];
            end
        end
    end
    %----------------------------------------------------------------------
    function mergeGroupsWithOnly1Difference()
        % If two electrodes have MORE THAN ONE electrode AND only one
        % electrode is different between them, merge.
        for i=length(groups):-1:1
            g1idx = groups{i};
            if length(g1idx) == 1
                continue
            end
            for k=i-1:-1:1
                g2idx = groups{k};
                if g2idx == 1
                    continue
                end
                a = setdiff(g1idx, g2idx);
                b = setdiff(g2idx, g1idx);
                if length(a)<=1 && length(b)<=1
                    D = max(computeDists([x(g1idx) y(g1idx)], [x(g2idx) y(g2idx)]),[],2);
                    if min(D) < P.addIfNearerThan
                        nGroupsPerElectrode(g1idx) = nGroupsPerElectrode(g1idx) -1;
                        nGroupsPerElectrode(g2idx) = nGroupsPerElectrode(g2idx) -1;
                        groups{k} = union(g1idx, g2idx);
                        nGroupsPerElectrode(groups{k}) = nGroupsPerElectrode(groups{k}) +1;
                        groups(i) = [];
                        break
                    end
                end
            end
        end
    end
    %----------------------------------------------------------------------
    function d = computeDists(p, X)
        nP = size(p,1);
        nX = size(X,1);
        d = zeros(nX,nP);
        for i=1:nX
            for k=1:nP
                d(i,k) = norm(p(k,:)-X(i,:));
            end
        end        
    end
end
