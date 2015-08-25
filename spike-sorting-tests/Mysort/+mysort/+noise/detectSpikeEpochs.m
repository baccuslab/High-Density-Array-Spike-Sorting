function spikeEpochs = detectSpikeEpochs(X, thresholds, minLen)
    X = mysort.datasource.DataObject(X);
    nC = size(X,1);
    spikeEpochs = cell(nC,1);
    if isscalar(thresholds)
        thresholds = ones(thresholds,1)*thresholds;
    end
    if nargin < 3
        minLen = 30;
    end
    
    for i=1:nC
       spikeEpochs{i} = mysort.epoch.fromBinaryVectorMinLen(abs(X(i,:))>=thresholds(i), minLen);
    end