function [clustCent, data2cluster, cluster2dataCell] = MeanShiftCluster2(X, bandWidth, varargin)
    % Mean-Shift Clustering of data using a flat kernel
    P.nMinPointsInCluster =[];
    P.maxNBandwidthIncreases = [];
    P.bandwidthIncreaseFactor = [];
    P.maxNGradientSpikes = 2 *1e6;
    P = mysort.util.parseInputs(P, varargin, 'error');
    

    % INIT
    [numPts, numDim] = size(X);
    numClust         = 0;
    bandSq           = bandWidth^2;
    pointEquivalence  = 1e-3*bandWidth;   % when points are closer than this, they are merged
    ids = (1:numPts)';
    
    
    % Make a copy of the points. This copy fill stay fixed and is used to
    % estimate the gradient at each point that is moving
    if P.maxNGradientSpikes < numPts
        rIdx = randperm(numPts);
        G = X(rIdx(P.maxNGradientSpikes),:);
    else
        G = X;
    end
    
    % ACTUAL ITERATION
    while 1
        X = ms(X);
    end
    
    function M = ms(X)
        
    end
end
