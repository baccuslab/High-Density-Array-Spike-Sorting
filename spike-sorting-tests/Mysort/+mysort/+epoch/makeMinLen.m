
function epochs = makeMinLen(epochs, minLen)
    if isempty(epochs)
        return
    end
    lens = epochs(:,2)-epochs(:,1)+1;
    diffs = minLen-lens;
    diffs(diffs<0) = 0;
    
    epochs(:,1) = epochs(:,1)-ceil(diffs/2);
    epochs(:,2) = epochs(:,2)+floor(diffs/2);