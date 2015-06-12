function epochs = fromBinaryVectorMinLen(binVec, minLen)
    epochs = [];
    if isempty(binVec)
        return
    end
    epochs = mysort.epoch.fromBinaryVector(binVec); 
    L = length(binVec);
    if isempty(epochs)
        return
    end
    epochs = mysort.epoch.makeMinLen(epochs, minLen);
    epochs = mysort.epoch.merge(epochs);
    epochs(epochs(:,1)<1,1) = 1;
    epochs(epochs(:,2)>L,2) = L;