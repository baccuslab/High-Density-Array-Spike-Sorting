function [FR, x] = toFiringRate(ST, window, sampling, timePeriod)
    % converts spike trains to vectors of firing rates
    % Input:
    %  ST       - cell array of spike Trains in SECONDS !
    %  window   - width of kernel for smoothing in seconds
    %  sampling - sampling (delta t) of the returned firing rate in
    %             seconds!
    %  timePeriod  - vector with two elements giving the temporal length of
    %                the resulting firing rate vectors as start and stop
    %                time (in seconds)
    %
    % Output:
    %  FR       - time series of firing rates if ST is a vector
    %             cell array of time series if ST is a cell array
    
    x = [];
%     if isempty(ST)
%         FR = ST;
%         return
%     end
    
    if iscell(ST)
        FR = cell(size(ST));
        for i = 1:size(ST,1)
            for j = 1:size(ST,2)
                [FR{i,j}, x] = mysort.spiketrain.toFiringRate(ST{i,j}, window, sampling, timePeriod);
            end
        end
        return
    else
        timePeriod = timePeriod(:);
        assert(length(timePeriod) == 2, 'must be 2 element vector');
        minTime = timePeriod(1);
        maxTime = timePeriod(2);
        nBins = ceil((maxTime-minTime)/sampling) + 1;
        x = (0:nBins-1)*sampling + minTime;
        FR = histc(ST, x);
        windowInSamples = ceil(window/sampling);
        tauRange = -2*windowInSamples:2*windowInSamples;
        kernel = normpdf(tauRange, 0, windowInSamples/2);
        kernelI = trapz(tauRange*sampling, kernel); % normalization?
        kernel = kernel/kernelI;
        FR = conv(FR, kernel, 'same');
    end
    
    
    
    
    
    