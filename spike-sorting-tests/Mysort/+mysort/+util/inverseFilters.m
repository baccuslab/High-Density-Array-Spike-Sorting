
function F = inverseFilters(XI,C)
    % XI: Tensor containing the templates (Tf x nC x nF)
    % C: Noise Covariancematrix
    
    minCond = 100;
    
    
    nC = size(XI,2);
    C = mysort.util.diagonalSubspaceLoading(C, minCond);
    iC = inv(C);
    XIconc = mysort.wf.t2m(XI);
    Fconc = XIconc*iC;
    F = mysort.wf.m2t(Fconc,nC); 