function [X tau] = vAlignOnTemplateInterpolated(X,nC,t,maxShift)
    if nargin < 4
        maxShift = 5;
    end
    opt = [];
    opt.MaxFunEvals = 30;
    opt.TolX = 0.01;
    opt.Display = 'off';    
    
    nS = size(X,1);

    mt = mysort.wf.v2m(t, nC);
    range = 1:size(mt,2);
    
    tau = zeros(nS, 1);
    for i=1:nS
        xsinc = mysort.wf.mSincfun(mysort.wf.v2m(X(i,:), nC));            
        fun = @(t)   -sum(sum(mt.*xsinc(range-t)));
        [tau(i), k] = fminbnd(fun, -maxShift, maxShift, opt);
    end

    X   = mysort.util.shiftRowsInterpolated(X, tau, nC);
end