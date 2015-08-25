function [aliV tau mx] = vAlignOnMax(V, nC, varargin)
    P.debug = false;
    P.truncate = 1;
    P.maxIdx = [];
    P.restrictToIdx = [];
    P = mysort.util.parseInputs(P, 'alignWaveforms', varargin);
    Tf = size(V,2)/nC;
    assert(round(Tf)==Tf, 'Number of channels does not match!');
    
    % convert to tensor
    T = mysort.wf.v2t(V, nC);
    % Check if we have to ignore some parts of the waveforms
    offset = 0;
    if ~isempty(P.restrictToIdx)
        idx = setdiff(1:size(T,1), P.restrictToIdx);
        T(idx,:,:) = 0;
        offset = idx(1)-1;
    end
    % compute max over dim 2 (channels)
    MT = squeeze(max(T,[],2));
    % find maximum over dim 1 (time points)
    [mx idx] = max(MT);
    idx = idx + offset;
    % mx contains the maxima of every waveform over all channels
    % idx contains the timepoint when the maxima occured. The channel
    % information is gone.
    
    if isempty(P.maxIdx)
        P.maxIdx = floor(size(MT,1)/2);
        tau = (-(idx - P.maxIdx))';
        % if there is an even number of wfs we dont want to end up with
        % the median in the middle of two values
        med = round(median(tau));
        tau = tau - med;   
    else
        tau = (-(idx - P.maxIdx))';
    end
    

    aliV = mysort.wf.vShift(V, nC, tau, P.truncate);
    
    if P.debug
        mysort.plot.spikes(V,'nC',nC);
        title('Raw Spikes');
        mysort.plot.spikes(MT,'nC',1);
        title('Preprocessed Spikes');            
    end
end