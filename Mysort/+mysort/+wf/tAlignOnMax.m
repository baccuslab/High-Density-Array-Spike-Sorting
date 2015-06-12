function [A tau mx] = tAlignOnMax(T, varargin)
    P.debug = false;
    P.truncate = 1;
    P.maxIdx = [];
    P.restrictToIdx = [];
    P.restrictToChannels = [];
    P = mysort.util.parseInputs(P, varargin, 'error');
    
    A = T;
    
    % Check if we have to ignore some parts of the waveforms
    offset = 0;
    if ~isempty(P.restrictToIdx)
        idx = setdiff(1:size(T,1), P.restrictToIdx);
        A(idx,:,:) = 0;
        offset = idx(1)-1;
    end
    if ~isempty(P.restrictToChannels)
        idx = setdiff(1:size(T,2), P.restrictToChannels);
        A(:,idx,:) = 0;
    end
    % compute max over dim 2 (channels)
    MT = squeeze(max(A,[],2));
    % find maximum over dim 1 (time points)
    [mx idx] = max(MT);
    idx = idx + offset;
    % mx contains the maxima of every waveform over all channels
    % idx contains the timepoint when the maxima occured. The channel
    % information is gone.
    
    if isempty(P.maxIdx)
        % if there is an even number of wfs we dont want to end up with
        % the median in the middle of two values
        med = round(median(idx));
        tau = -(idx - med);
    else
        tau = (-(idx - P.maxIdx))';
    end
    tau = tau(:);
    A = mysort.wf.tShift(T, tau, P.truncate);
end