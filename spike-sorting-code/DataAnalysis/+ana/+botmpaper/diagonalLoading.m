function [C] = diagonalLoading(C, condNumber)
    % inputs:
    %    C - Covariance matrix
    %    condNumber - minimal condition number, C should have
    % 
    % algorithm:
    % It should be EVmax/EVmin = condNumber
    
    [V, D] = eig(C);
    D = diag(D);
    maxEV = max(D);
    minEV = min(D);
    assert(condNumber>=1, 'Condition number cannot be smaller than 1!');
    
    if maxEV == 0 || condNumber == 1
        C = eye(size(C));
        return
    end
    
    assert(minEV>0, 'This should only be applied to semi positive definit matrices!');
    
    if maxEV/minEV <= condNumber
        return
    end
        
    lambda = (maxEV-condNumber*minEV)/(condNumber-1);
    
    C = C + lambda*eye(size(C));
end    