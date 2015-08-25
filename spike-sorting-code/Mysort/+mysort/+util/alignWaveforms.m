
function [nettoShift aliX] = alignWaveforms(X, nC, varargin)
    warning('This function is depricated! Use mysort.wf.* instead!');
    P.nIter = 3;
    P.debug = false;
    P.truncate = 1;
    P = mysort.util.parseInputs(P, varargin);
%     mysort.plot.spikes(X,'nC',nC); title('BEFORE');
    [MX nettoShift] = mysort.wf.vAlignOnMax(X, nC, 'truncate', 1);
%     mysort.plot.spikes(MX,'nC',nC); title('AFTER');

    for iter=1:P.nIter
        tau = alignStep(MX, nC);
        MX = mysort.util.shiftMCRows(MX,tau,nC,0);
        nettoShift = nettoShift + tau;  
        nShifts = sum(tau~=0);
        fprintf('alignWaveforms: %d shifted\n', nShifts);
        if nShifts == 0
            break;
        end
    end
    
    if nargout > 1
        aliX = mysort.util.shiftMCRows(X, nettoShift, nC, P.truncate);
    end
    
    function tau_ = alignStep(X, nC)
        mx = mean(X,1);
        yx = zeros(size(X));
        for i=1:size(X,1)
            yx(i,:) = mysort.util.mcfilt(X(i,:), mx);
        end
        [m tau_] = max(yx,[],2);
        tau_ = -(tau_ - ceil(length(mx)/2));
        tau_ = tau_ - median(tau_);

        if P.debug
            mysort.plot.spikes(X,'nC',nC);
            title('Raw Spikes');
            mysort.plot.spikes(XX,'nC',1);
            title('Preprocessed Spikes');
            mysort.plot.spikes(mx,'nC',1);
            title('Mean = Filter');
            mysort.plot.spikes(yx,'nC',1);
            title('Filter Outputs');            
        end
    end
end