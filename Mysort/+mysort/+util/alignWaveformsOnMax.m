
function [tau aliX] = alignWaveformsOnMax(X, nC, varargin)
warning('This function is depricated! Use mysort.wf.* instead!');
    P.debug = false;
    P.truncate = 1;
    P.maxIdx = [];
    P.restrictToIdx = [];
    P = mysort.util.parseInputs(P, 'alignWaveforms', varargin);
    Tf = size(X,2)/nC;
    assert(round(Tf)==Tf, 'Number of channels does not match!');
    
    XX = mysort.wf.v2t(X, nC);
    if ~isempty(P.restrictToIdx)
        idx = setdiff(1:size(XX,1), P.restrictToIdx);
        XX(idx,:,:) = 0;
    end
    MX = squeeze(max(XX,[],2));
    [mx idx] = max(MX);
    
    if isempty(P.maxIdx)
        P.maxIdx = floor(size(MX,1)/2);
        tau = (-(idx - P.maxIdx))';
        tau = tau - median(tau);        
    else
        tau = (-(idx - P.maxIdx))';
    end
    

    if nargout > 1
        aliX = mysort.util.shiftMCRows(X, tau, nC, P.truncate);
    end
    if P.debug
        mysort.plot.spikes(X,'nC',nC);
        title('Raw Spikes');
        mysort.plot.spikes(MX,'nC',1);
        title('Preprocessed Spikes');            
    end
end