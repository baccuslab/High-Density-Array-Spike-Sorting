function GT = E10_preprocessing(load_path, gtFile, cutLeft, Tf, save_path)
% this is so fast, dont buffer.
%     bufferFile = fullfile(save_path, gtFile, 'preproc.mat');
%     if exist(bufferFile, 'file')
%         fprintf('Loading from buffer ... %s\n', bufferFile);
%         load(bufferFile)
%         return
%     end
    fprintf('Preprocessing ... %s\n', gtFile);
    data = load(fullfile(load_path, gtFile));
    upsample = 1;    % Dont upsample the whole data. That is far too slow.
    cutLeft = cutLeft * upsample;
    Tf = Tf * upsample;
    
    % The struct GT will contain all the necessary information for the 
    % evaluation and is returned
    GT = [];
    GT.source_filename = gtFile;
    GT.X = data.data; % upsampling far too slow!!! mysort.util.resampleMC(data.data, upsample, 1);
    GT.srate = 24000*upsample;
    GT.nC = size(GT.X,1);
    GT.gdf(:,1) = data.spike_class{1};
    GT.gdf(:,2) = data.spike_times{1}-cutLeft;
    GT.classes = unique(GT.gdf(:,1));
    GT.spikeEpochs = [data.spike_times{1}'-cutLeft data.spike_times{1}'-cutLeft+Tf-1];
    GT.spikeEpochs = GT.spikeEpochs(GT.spikeEpochs(:,1)>0,:);
    GT.spikeEpochs = GT.spikeEpochs(GT.spikeEpochs(:,2)<size(GT.X,2),:);  
    
    % extract Spikes   
    GT.spikesX = mysort.epoch.extractWaveform(GT.X, GT.spikeEpochs);
    
    [GT.noiseEpochs GT.C] = mysort.util.calculateNoiseEpochsAndCovarianceMatrix(...
                                    GT.X, GT.spikeEpochs, Tf);
    XI = [];
    for f=1:length(GT.classes)
        XI(:,1,f) = mean(GT.spikesX(GT.gdf(:,1)==GT.classes(f),:));
    end
    GT.XIunaligned = XI;
    GT.tau = zeros(1, size(XI,3));
    [absXI GT.tau] = mysort.wf.tAlignOnMax(abs(XI), 'truncate', 1); 
    XI = mysort.wf.tShift(XI, GT.tau, 1);
%     [XI GT.tau] = mysort.wf.tAlignOnCorrelation(XI, 'trunc', 1);
    
    GT.XI = XI;
    GT.templates = mysort.wf.t2m(XI);
    fprintf('...preprocessing done.\n');
    