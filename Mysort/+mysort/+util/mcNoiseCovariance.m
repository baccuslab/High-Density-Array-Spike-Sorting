
function C = mcNoiseCovariance(X, noiseEpochs, Tf, minEpochLength)
    if nargin == 2
        % If only 2 parameters are supplied, X is treated as only noise
        % and only Tf has to be supplied.
        Tf = noiseEpochs;
        noiseEpochs = [1 size(X,2)];
    end

    % if minimal noise epoch length is not provided use default of 10.
    if nargin < 4
        minEpochLength = 10;
    end
    
    nC = size(X,1);
    C = zeros(nC*Tf);
    
    epochLength = noiseEpochs(:,2)-noiseEpochs(:,1)+1;
    noiseEpochs = noiseEpochs(epochLength>minEpochLength,:);
    epochLength = epochLength(epochLength>minEpochLength);
    
    totalNoiseEpochLength = sum(epochLength);
    for i=1:size(noiseEpochs,1)
        C = C + epochLength(i,1) * mysort.util.mcDataCovariance(X(:, noiseEpochs(i,1):noiseEpochs(i,2)), Tf);
    end
    C = C/totalNoiseEpochLength;