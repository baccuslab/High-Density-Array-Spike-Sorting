function fp = findPointsInEpochs(points, epochs)
    % finds all points in points that happen during an epoch in epochs
    pointPtr = 1;
    epochPtr = 1;
    
    points = sort(points);
    epochs = sortrows(epochs);
    
    nP = length(points);
    nE = size(epochs,1);
    
    fp = zeros(nP,1);
    while pointPtr<=nP && epochPtr<=nE
        if points(pointPtr) < epochs(epochPtr,1)
            % current point is before epoch
            pointPtr = pointPtr+1;
        elseif points(pointPtr) <= epochs(epochPtr,2)
            % current point is in current epoch
            fp(pointPtr) = epochPtr;
            pointPtr = pointPtr+1;
        else
            % current point is behind current epoch
            epochPtr = epochPtr+1;
        end        
    end        
end