% function [sorting SO] = meanSpike(X, Tf)
fprintf('Init...\n');
    [nC Len] = size(X);
    thr_noise = 4;
    thr_spike = 4.5;
    cutleft = floor(Tf/3);
    cutright = Tf-1-cutleft;
    
%%    
fprintf('Estimate channel wise noise std and cut single channel spikes...\n');    
    [smad, noiseepochs] = mysort.noise.estimateSigma(X, Tf, thr_noise);
    SO.init.noiseepochs = noiseepochs;
    SO.init.singleChannelSpikeTimes = {};
    SO.init.singleChannelSpikes = {};
    SO.init.singleChannelSpikeEpochs = {};
%%
    X = X./repmat(smad, 1, L); 
%%    
fprintf('Cutting spike waveforms on every channel...\n');  
    chunker = mysort.util.Chunker(nC, 'chunkSize', 1, 'progressDisplay', 'console');
    while chunker.hasNextChunk()
        i = chunker.getNextChunk(); i = i(1);
        [pks, times] = findpeaks(X(i,:), 'MINPEAKHEIGHT', thr_spike,...
                                  'MINPEAKDISTANCE', Tf);
        epochs = [times'-cutleft times'+cutright];
        SO.init.singleChannelSpikeTimes{i} = times;
        SO.init.singleChannelSpikeEpochs{i} = epochs;
        SO.init.singleChannelSpikes{i} = mysort.epoch.extractWaveform(X(i,:),epochs);
    end

    %% ALIGN SPIKES
fprintf('Align single channel spikes and calc mean spike...\n');
    scs = SO.init.singleChannelSpikes;
    allSpikes = cell2mat(scs');
    [nettoShift alignedSpikes] = mysort.util.alignWaveformsInterpolatedMean(allSpikes, 1);

    %% CONSTRUCT big Spike mean
    mX = mean(alignedSpikes, 1);
    mX = mX/norm(mX);
    
    %% group spike times
fprintf('Group single channel spikes to multichannel spikes...\n');    
    SO.init.spikeTimes = sort([SO.init.singleChannelSpikeTimes{:}]);
    lastspike = -100;
    remove = zeros(length(SO.init.spikeTimes), 1);
    for i=1:length(SO.init.spikeTimes)
        if SO.init.spikeTimes(i)-lastspike < Tf
            remove(i) = 1;
        else
            lastspike = SO.init.spikeTimes(i);
        end
    end
    SO.init.spikeTimes(remove==1) = [];

    %% EXTRACT FULL WAVEFORMS
fprintf('Extract multichannel spikes...\n');    
    SO.init.spikeEpochs = [SO.init.spikeTimes'-cutleft SO.init.spikeTimes'+cutright];
    spikes = mysort.epoch.extractWaveform(X, SO.init.spikeEpochs);

    %% PROJECTION ON MEAN
fprintf('Project on meanspike...\n');    
    XX = mysort.util.v2m(spikes, nC);
    XXm = max(mysort.util.convMatrix(XX, mX), [], 2);
    %XXm = XX*mX';
    Xm = mysort.util.m2v(XXm, nC);

    %% remove "empty" channels
    maxi = max(Xm);
    Xm_red = Xm(:,maxi>200);
    if 0
        figure; plot(Xm_red', '.')
        figure; plot(Xm'); hold on; plot(maxi, 'xm')    
        figure;
        ah(1) = subplot(1,2,1);
        plot3(Xm_red(:,16), Xm_red(:,34), Xm_red(:,46), '.')
        linkdata on
        ah(2) = subplot(1,2,2);
        plot(Xm_red', '.')
        linkdata on
%         databrushing on
        
        
    end
    %% PCA 
fprintf('PCA feature extraction...\n');    
    FET = mysort.util.dimReductionPCA(Xm_red, 6);
    if 0
        figure; plot(FET', '.')
        
        figure; plot3(FET(:,1), FET(:,2), FET(:,3), '.')
    end    
    %% CLUSTERING
fprintf('K Means clustering...\n');    
%     [classes centers N obj R] = mysort.clustering.gmmFit(X, ...
%     'kmin', 1, 'kmax', 40, 'prewhitened', false, 'homoscedastic', true);
% [classes BIC minK] = mysort.clustering.kmeansBIC(FET, 'kmin', 8, 'kmax', 25, 'repeats', 50);
[classes centers N obj R] = mysort.clustering.gmmClustering(FET, ...
    'repeats', 50, ...
    'kmin', 1, 'kmax', 35, 'prewhitened', true, 'homoscedastic', true,...
    'fixedCov', ones(1, size(FET,2))); %[]);

    %% COMPUTE TEMPLATES
fprintf('Compute templates...\n');    
    cl = unique(classes);
    T = zeros(length(cl), nC*Tf);
    n = zeros(length(cl),1);
    for i=1:length(cl)
        T(i,:) = mean(spikes(classes==cl(i),:),1);
        n(i)   = sum(classes==cl(i));
    end

    %% OUTPUT
fprintf('Prepare output...\n');    
    SO.sorting = [classes SO.init.spikeTimes'];
    if nargout > 1
        SO.templates = T;
        SO.clustering.data = Xm;
        SO.clustering.features = FET;
        SO.meanSpike = mX;
        SO.spikes = spikes;
        SO.clustering.BIC = BIC;
        SO.clustering.minK = minK;
        SO.init.thr_noise = thr_noise;
        SO.init.thr_spike = thr_spike;
        SO.init.smad = smad;
        SO.init.cutleft = cutleft;
        SO.init.cutright = cutright;
    end
    sorting = SO.sorting;
