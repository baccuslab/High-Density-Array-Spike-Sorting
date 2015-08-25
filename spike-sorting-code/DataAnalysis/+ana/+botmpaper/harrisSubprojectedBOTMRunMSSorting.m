def = mysort.util.defs();
names = [1:4 6];
for nameidx = 5:length(names)
    name = names(nameidx)
     fprintf('Starting with %d\n', name);
   close all
   pic = 0;
    if 1
    %     name = 2;
        [D X Xfil I Ifil] = ana.harris.preprocess.preprocess(name);
        [pks gtSpikeTimes] = findpeaks(I+.0001*randn(size(I)), 'minpeakheight', -2000, 'minpeakdistance', 30);
    
        if 0
            figure; plot(I);
            hold on
            plot(gtSpikeTimes, pks, 'rx', 'markersize', 20, 'linewidth', 3);
        end
        
        DS = mysort.ds.Matrix(Xfil', D.srate, 'bla');
        dpath = fullfile(D.path, '..', 'MeanShiftSorting6', D.name);

        P.spikeDetection.method = '-';
        P.spikeDetection.thr = 4.0;
        P.spikeDetection.minDist = 35;
        
        P.spikeCutting.Tf = 55;
        P.spikeCutting.cutLeft = 15;
        P.spikeCutting.maxSpikes = 200000;
        
        % Spike Alignment
        P.spikeAlignment.method = 'onUpsampledMean';
        P.spikeAlignment.maxSpikes = 200000;
        P.spikeAlignment.Tf = 25;
        P.spikeAlignment.cutLeft = 10;
        P.spikeAlignment.initAlignment = '-';
        P.spikeAlignment.maxIdx = P.spikeAlignment.cutLeft + 1;
        P.spikeAlignment.maxIterations = 30;

        % Noise estimation    
        P.noiseEstimation.minLength = 2000000;
        P.noiseEstimation.minDistFromSpikes = 60;

        % Feature extraction
        P.featureExtraction.Tf = 14;  
        P.featureExtraction.cutLeft = 7;
        P.featureExtraction.nDims = 6;
    
        % Clustering
        P.clustering.minSpikesPerCluster = 10;
        P.clustering.meanShiftBandWidth = sqrt(1.2*P.featureExtraction.nDims);
        P.clustering.maxSpikes = P.spikeAlignment.maxSpikes;
        
        P.botm.cutLeft = P.spikeCutting.cutLeft;
        
        [S P] = ana.sort(DS, dpath, 'r2', P)
            
        Tf = P.spikeCutting.Tf;
    end
end