function [Xstar radVecs S]= starTransform(X, d)
    if nargin < 2
        d = 2;
    elseif d ~= 2
        error('not implemented');
    end
    
    if d==2
        nS = size(X,d);
        radVecs = mysort.util.getRadialVectors(nS, d);
%         % Find weights of inter radial vector similarity (aka scalar
%         % product)
        W = abs(radVecs* radVecs');
        % Compute covariance and ignore diagonal
        C = cov(X) + diag(nan(1, nS));      
        % Find best permutation of C so that similar (according to C and
        % weighted by W) columns are next to each other (circular)
        S = mysort.util.findBestWeightedPermutation(C, W);
        %S = 1:nS;
        % project data accoring to found order
        Xstar = zeros(size(X,1),d);
        for i=1:nS
            k = S(i);
            Xstar = Xstar + X(:,k)*radVecs(i,:);
        end
    end
    
end