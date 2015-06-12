
function [tau aliX] = vAlignOnAverageMaxSample(X, nC, varargin)
    P.debug = false;
    P.truncate = 1;
    P.maxIdx = [];
    P.restrictToIdx = [];
    P.nMaxChannelsForWeighting = 5;
    P = mysort.util.parseInputs(P, 'alignWaveforms', varargin);
    Tf = size(X,2)/nC;
    assert(round(Tf)==Tf, 'Number of channels does not match!');
    
    XX = mysort.wf.v2t(X, nC);
    if ~isempty(P.restrictToIdx)
        idx = setdiff(1:size(XX,1), P.restrictToIdx);
        XX(idx,:,:) = 0;
    end
    
    %Compute for each waveform the position of the maximum on every channel
    [MX maxSamplePerChannel] = max(XX,[],1);
    MX = squeeze(MX);
    maxSamplePerChannel = squeeze(maxSamplePerChannel);
    
    % weight with the amplitudes 
    % if one channel has less then 50% amp than the max channel, ignore
    allChannelMax = max(MX,[],1);
    WMX = MX;
    for i=1:nC
        WMX(i,MX(i,:)<allChannelMax*.5) = 0;
    end
    
    % allow only the nMaxChannelsForWeighting biggest channels for the
    % weighting
    if P.nMaxChannelsForWeighting < nC
        for i=1:size(WMX,2)
            spikeHeights = sortrows([WMX(:,i) (1:nC)']);
            WMX(spikeHeights(1:(end-P.nMaxChannelsForWeighting),2),i) = 0;
        end
    end
    
    % normalize to 1
    WMX = WMX./ repmat(sum(WMX,2), nC, 1);
    weightedMaxSamplePerChannel =  sum(WMX.*maxSamplePerChannel,2);
    idx = weightedMaxSamplePerChannel;
    if isempty(P.maxIdx)
        P.maxIdx = floor(size(MX,2)/2);
        tau = (-(idx - P.maxIdx))';
    else
        tau = (-(idx - P.maxIdx))';
    end
    
    if nargout > 1
%         aliX = mysort.util.shiftMCRows(X, round(tau), nC);
        XX = mysort.wf.v2t(X, nC);
        aliX = mysort.wf.t2v(XX);
        aliX = mysort.util.shiftRowsInterpolated(aliX, tau , nC);
    end
    if P.debug
        mysort.plot.spikes(X,'nC',nC);
        title('Raw Spikes');
        mysort.plot.spikes(MX,'nC',1);
        title('Preprocessed Spikes');            
    end
end