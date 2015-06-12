function [std, mad] = mad(X)
% Calculates the MAD (median absolute deviation) for every row of X.
    k = 1.4826;
    nC = size(X,1);
    
    mad = zeros(nC,1);
    for c = 1:nC
        med = median(X(c,:));
        mad(c) = median(abs(X(c,:)-med));
    end
    std = mad*k;