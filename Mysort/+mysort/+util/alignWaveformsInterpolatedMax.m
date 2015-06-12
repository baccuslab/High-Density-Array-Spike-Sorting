
function [tau S] = alignWaveformsInterpolatedMax(S,nC, varargin)
    warning('This function is depricated! Use mysort.wf.* instead!');
    P.maxIdx = [];
    P = mysort.util.parseInputs(P, 'alignWaveformsInterpolatedMean', varargin);
    
    if nargin < 2 
        nC = 1;
    end
    
    Tf = size(S,2)/nC;

    opt = [];
    opt.MaxFunEvals = 25;
    opt.TolX = 0.05;
    opt.Display = 'off';    

    fprintf('mean variance: %f\n', mean(var(S)));
    maxpos = alignStep(S, nC);
    if isempty(P.maxIdx)
        P.maxIdx = median(maxpos);
    end
    tau = P.maxIdx-maxpos;    
    
    S   = mysort.util.shiftRowsInterpolated(S, tau, nC);
    nShifts = sum(tau~=0);
    fprintf('alignWaveforms: %d shifted\n', nShifts);
    fprintf('mean variance: %f\n', mean(var(S)));
    
    function maxpos = alignStep(X, nC)

        maxpos = zeros(size(X,1), 1);
        for i=1:size(X,1)

            xsinc = mysort.wf.mSincfun(mysort.wf.v2m(X(i,:), nC));            
            fun = @(t) -max(max(xsinc(t)));
            [m maxi] = max(X(i,:));
            maxi = mod(maxi-1, Tf)+1;
            %[tau(i), k] = fminsearch(fun, 0, opt);
            [maxpos(i), k] = fminbnd(fun, max(1,maxi-1), min(Tf, maxi+1), opt);

        end
    end
end