
function e = intersect(e1, e2, dont_sort, dont_merge)
    % Finds the common epochs in e1 and e2
    % if you are sure e1 and e2 are already sorted, provide dont_sort
    % if you are sure e1 and e2 are already merged, i.e. epochs within
    % one of the sets do not overlap, provide dont_merge
    if nargin < 3
        e1 = sortrows(e1);
        e2 = sortrows(e2);
    end
    if nargin < 4
        e1 = mysort.epoch.merge(e1);
        e2 = mysort.epoch.merge(e2);
    end
    e = [];
    if isempty(e1) || isempty(e2)
        return
    end
    i=1; k=1; 
    count = 1;
    e = zeros(max(size(e1,1),size(e2,1)), 2);
    while i<=size(e1,1) && k<=size(e2,1)
        % check if i-th epoch of e1 overlaps with k-th epoch of e2
        
        % Check if e1 is too far left, then iterate i
        % e1:  +++++
        % e2:         ++++++ 
        if e1(i,2) < e2(k,1)
            i = i+1;
        
        % Check if e2 is too far left, then iterate k
        % e1:         +++++
        % e2:  ++++++ 
        elseif e1(i,1) > e2(k,2)
            k=k+1;
            
        % Epochs are overlapping with one of these cases:            
            %         A           B        C      D
            % e1:  +++++++     +++++++    +++      +++++        
            % e2:     ++++++     ++++    +++++   ++++  
        else
            e(count,:) = [max(e1(i,1), e2(k,1)) min(e1(i,2), e2(k,2))];
            if e1(i,2) <= e2(k,2)
                i = i+1;
            else
                k = k+1;
            end                   
            count = count +1;
        end
    end
    e = e(1:count-1,:);