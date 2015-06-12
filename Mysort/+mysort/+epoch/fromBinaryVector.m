
function epochs = fromBinaryVector(binVec) 
    if size(binVec,2) == 1
        binVec = binVec';
    end
    starts = [binVec(1) binVec(2:end)   & ~binVec(1:end-1)];
    ends   = [binVec(1:end-1) & ~binVec(2:end) binVec(end)];

    epochs = [find(starts(:)) find(ends(:))];
end             