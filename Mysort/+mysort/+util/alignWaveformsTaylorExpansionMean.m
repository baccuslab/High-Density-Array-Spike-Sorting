
function [nettoShift S] = alignWaveformsTaylorExpansionMean(S,nC, varargin)
warning('This function is depricated! Use mysort.wf.* instead!');
    P.nIter = 3;
    P = mysort.util.parseInputs(P, 'alignWaveformsTaylorExpansionMean', varargin);
    
    if nargin < 2 
        nC = 1;
    end

    opt = [];
    opt.MaxFunEvals = 30;
    opt.TolX = 0.01;
    opt.Display = 'off';    
    nettoShift = zeros(size(S,1), 1);
    last_tau  = ones(size(S,1), 1);
    fprintf('mean variance: %f\n', mean(var(S)));
    for iter=1:P.nIter
        [tau last_tau] = alignStep(S, nC, last_tau);
        
        S   = mysort.util.shiftRowsInterpolated(S, tau, nC);
        nettoShift = nettoShift + tau;  
        nShifts = sum(tau~=0);
        fprintf('alignWaveforms: %d shifted\n', nShifts);
        fprintf('mean variance: %f\n', mean(var(S)));
        if nShifts == 0
            break;
        end
    end
    
    function [tau last_tau] = alignStep(X, nC, last_tau)
        mx = mysort.util.v2m(mean(X,1), nC);
        range = 1:size(mx,2);
        tau = zeros(size(X,1), 1);
        for i=1:size(X,1)
            if abs(last_tau(i))>.3
                xsinc = mysort.util.sincfun(mysort.util.v2m(X(i,:), nC));            
                fun = @(t)   -sum(sum(mx.*xsinc(range-t,0)));
                %[tau(i), k] = fminsearch(fun, 0, opt);
                [tau(i), k] = fminbnd(fun, -5, 5, opt);
                last_tau(i) = tau(i);
            else
                last_tau(i) = 1;
            end
%             fprintf('%d', i);
        end
%         fprintf('\n');
%         tau = tau-mean(tau);
    end
end