function [times pks] = detecSingleElectrodeSpikes(dataSource, varargin)
    if isnumeric(dataSource)
        dataSource = mysort.ds.Matrix(dataSource);
    end
    P.channelIdx = 1:dataSource.size(2);
    P.chunkSize = 100000;
    P.thr = 3.5;
    P.noiseStd = [];
    P.minPeakDistance = [];
    P = mysort.util.parseInputs(P, varargin);
    Len = dataSource.size(1);
    nC = length(P.channelIdx);
    if isempty(P.minPeakDistance)
        sps = dataSource.getSamplesPerSecond();
        if isempty(sps)
            P.minPeakDistance = 20;
        else
            P.minPeakDistance = ceil(sps/1000); % one ms
        end
    end

    % get noise std
    if isempty(P.noiseStd)
        disp('Estimating noise std...');
        X = dataSource(1:300000, P.channelIdx);
        [smad, commonNoiseEpochs, notNoiseEpochsSet] = mysort.noise.estimateSigma(X, P.minPeakDistance*2, P.thr);
        P.noiseStd = smad;
    elseif length(P.noiseStd) == 1 && nC>1
        P.noiseStd = repmat(P.noiseStd, 1, nC);
    end
    assert(length(P.noiseStd) == nC, 'Number of noise stds does not match requested channel count!');
    

    % Detect spikes in the beginning of the file
    disp('Detecting spikes...'); tic
    pks = cell(length(P.channelIdx),1);
    times = pks;
    for cidx = 1:length(P.channelIdx)
        c = P.channelIdx(cidx);
        pks{c,1} = [];
        times{c,1} = [];
    end
    chunker = mysort.util.Chunker(Len, 'chunkSize', P.chunkSize, ...
        'progressDisplay', 'console', 'minChunkSize', 1000, 'chunkOverlap', 2*P.minPeakDistance);
    while chunker.hasNextChunk()
        [chunkOvp chunk] = chunker.getNextChunk();
        X = double(dataSource.getData(chunkOvp(1):chunkOvp(2), P.channelIdx));
        for cidx = 1:length(P.channelIdx)
            c = P.channelIdx(cidx);
            [pks_, times_] = findpeaks(abs(X(:,c)), 'MINPEAKHEIGHT', P.noiseStd(cidx)*P.thr,...
                                                    'MINPEAKDISTANCE', P.minPeakDistance);
            pks_ = X(times_,c); % get the right amplitudes! (sign!)
            pks_ = pks_(:);
            times_ = times_(:);
            % remove spikes that are outside this chunk
            rmvIdx = (times_+chunkOvp(1) < chunk(1)) | (times_+chunkOvp(1) > chunk(2));
            pks_(rmvIdx) = [];
            times_(rmvIdx) = [];

            pks{c,1} = [pks{c}; pks_];
            times{c,1} = [times{c}; times_+chunkOvp(1)-1];
        end
    end
    disp('Done.'); toc    
end