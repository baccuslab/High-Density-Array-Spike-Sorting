
function epochs = makeMaxLen(epochs, maxLen)
    if isempty(epochs)
        return
    end
    lens = epochs(:,2)-epochs(:,1)+1;
    diffs = lens-maxLen;
    diffs(diffs<0) = 0;
    
    epochs(:,1) = epochs(:,1)+ceil(diffs/2);
    epochs(:,2) = epochs(:,2)-floor(diffs/2);