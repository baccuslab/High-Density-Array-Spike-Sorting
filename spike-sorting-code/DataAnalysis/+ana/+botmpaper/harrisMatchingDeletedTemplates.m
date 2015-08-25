if ~exist('Xfil', 'var')
    name = 1;
    [D X Xfil I Ifil] = ana.harris.preprocess.preprocess(name);
    %%
    P = struct();
    P.spikeDetection.method = '-';
    P.spikeDetection.thr = 3.9;
    P.artefactDetection.use = 0;
    P.botm.run = 0;
    P.spikeCutting.maxSpikes = 200000000000; % Set this to basically inf

    P.noiseEstimation.minLength = 300000; % This also effects computation time
    P.noiseEstimation.minDistFromSpikes = 80;

    P.spikeAlignment.initAlignment = '-';   
    P.spikeAlignment.maxSpikes = 60000;     % so many spikes will be clustered
    P.clustering.maxSpikes = P.spikeAlignment.maxSpikes;  % dont align spikes you dont cluster...
    P.clustering.meanShiftBandWidth = sqrt(1.8*6);
    P.mergeTemplates.merge = 1;
    P.mergeTemplates.upsampleFactor = 3;
    P.mergeTemplates.atCorrelation = .93; % DONT SET THIS TOO LOW! USE OTHER ELECTRODES ON FULL FOOTPRINT TO MERGE 
    P.mergeTemplates.ifMaxRelDistSmallerPercent = 30;
    %%
    DS  = mysort.ds.Matrix(Xfil', D.srate, 'Prefiltered Extracellular');
    DSI = mysort.ds.Matrix(I', D.srate,    'Raw Intracellular');
    % dpath = ['C:\LocalData\HarrisBuzsaki\MeanShiftSorting\' D.folder];
    dpath = ana.harris.datapath();
    dpath = fullfile(dpath, '..', 'MeanShiftSorting8', D.name);
    [S P] = mysort.sorters.sort(DS, dpath, 'r1', P)
    Tf = P.spikeCutting.Tf;
    %%
    % mysort.plot.waveforms(S.clusteringMatched.templates)
    mysort.plot.waveforms(S.clusteringMatched.spikeCutAligned,...
        'IDs', S.clusteringMatched.ids, 'stacked', 0)
    
    sgdf = [S.clusteringMatched.ids(:) S.clusteringMatched.ts(:)];
    tgdf = [ones(length(D.gt_spiketrain),1) D.gt_spiketrain];
    R = mysort.spiketrain.alignGDFs(sgdf, tgdf, 10, 10, 10);
    mysort.plot.printEvaluationTable(R, 'ignoreEmptyGT', false, 'ignoreEmptySorted', false)

    %%
    defs = mysort.util.defs();
    doNotLabel = [defs.tp defs.tpo defs.fp defs.fn defs.fno];

    %%
    vT = D.gt_template;
    nSp = size(D.gt_spikes,1);
    Eucl = sum((D.gt_spikes-repmat(vT, nSp, 1)).^2,2);
    Conv = D.gt_spikes*vT';

    %% Noise Analysis
    spikeEpochs = mysort.epoch.merge([S.spikeDetectionMerged.ts(:)-80 S.spikeDetectionMerged.ts(:)+80]);
    noiseEpochs = mysort.epoch.flip(spikeEpochs, size(DS,1));
    noiseEpochs = mysort.epoch.removeShort(noiseEpochs, P.spikeCutting.Tf);
    noiseSnippets = DS.getWaveform(noiseEpochs(:,1), 0, P.spikeCutting.Tf);
    nNoise = size(noiseSnippets,1);
    
    T = S.clusteringMatched.templates;
    Tfbotm = size(T,2)/size(DS,2);
    LN = 1000000;
    Cest = mysort.noise.Covest2(DS, 'maxLag', Tfbotm-1, 'maxSamples', LN, 'noiseEpochs', noiseEpochs);
    C = mysort.noise.ccol2Cte(Cest.CCol);

    %%
    correctIdx = R.f2k;
    correctId = R.St1IDs(correctIdx);
    correctSpikes = DS.getWaveform(D.gt_spiketrain, P.spikeCutting.cutLeft+7, P.spikeCutting.Tf);
    Tcorrect = sum(mysort.wf.v2m(mean(correctSpikes), size(DS,2)));
    
    % Compute offset between ground truth and sorted spike train
    [mi miidx]  = min(Tcorrect);
    [mi miidx2] = min(sum(mysort.wf.v2m(T(correctIdx,:), size(DS,2))));
    offset = miidx-miidx2;
    
    % RECUT ALL SPIKES ACCORDING TO SORTING
    allSpikes = DS.getWaveform(round(sgdf(:,2)), P.spikeCutting.cutLeft+7-1, P.spikeCutting.Tf);
    
    % Compute Templates
    correctSpikes_gt = DS.getWaveform(D.gt_spiketrain, P.spikeCutting.cutLeft+7-offset, P.spikeCutting.Tf);
    Tcorrect_gt = mean(correctSpikes_gt, 1);
    
    U = unique(sgdf(:,1));
    correctSpikes_sorted_recut = {};
    Tcorrect_sorted_recut = [];
    for i=1:length(U)
        correctSpikes_sorted_recut{i} = DS.getWaveform(round(sgdf(sgdf(:,1)==U(i),2)), P.spikeCutting.cutLeft+7-1, P.spikeCutting.Tf);
        Tcorrect_sorted_recut(i,:) = mean(correctSpikes_sorted_recut{i});    
    end
        
end


%%
T_correct = Tcorrect_sorted_recut(correctIdx,:);
T_others  = Tcorrect_sorted_recut;
T_others(correctIdx,:) = [];
nOthers = size(T_others,1);
nRep = 100;
nDeleteTemplates = 0:7;
Sensitivity = [];
Specificity = [];
Classification = [];
T1 = [];
T2 = [];
for nd = 1:length(nDeleteTemplates)
    d = nDeleteTemplates(nd);
    for r = 1:nRep
        rp = randperm(nOthers);
        T_this = [T_correct; T_others(rp(1:nOthers-d),:)];
        RR = ana.botmpaper.harrisMatchingEvalFun(T_this, allSpikes, correctSpikes_sorted_recut{correctIdx}, noiseSnippets, C);
        correctIdx_this = 1;
        [mx mxidx] = max(RR.botmCor,[],2);
        nClassified = sum(mxidx==correctIdx_this & mx > 0);
        nCorrectSpikes = length(mxidx);
        classPerfBotm = nClassified/nCorrectSpikes;
        Classification(nd,r) = classPerfBotm;
        
        [mx mxidx] = max(RR.botm,[],2);
        % This is the matching to all spikes, so also the correctly matched
        % correct spikes. Thus, subtract those from this number 
        rejPerfBotm = (sum(mxidx==correctIdx_this & mx > 0)- nClassified)/(length(mxidx)-nCorrectSpikes);
        Specificity(nd,r) = rejPerfBotm;
        
        [mx mxidx] = max(RR.botmNoise,[],2);
        nClassified = sum(mxidx==correctIdx_this & mx > 0);
        nTotalNoise = length(mxidx);
        percFP = nClassified/nTotalNoise;
        Sensitivity(nd,r) = percFP;    
    end
end
nC = size(DS,2);
save('harrisMatchingDeletedTemplatesBuffer2.mat', 'Sensitivity', 'Classification', 'Specificity',...
    'nDeleteTemplates', 'correctIdx', 'nC');

return
%%
figure
plot(mean(allSpikes))
hold on
plot(mean(S.spikeCut.wfs), '.-g')
legend('Recut Spikes', 'Spikes From Sorter');
%%
figure
plot(Tcorrect_gt', '-bx')
hold on
plot(Tcorrect_sorted_recut(correctIdx,:)', 'g-d')
plot(T(correctIdx,:), 'c');
legend('Ground Truth', 'Sorted recut', 'Sorted from sorter');
%%
figure
plot(median(correctSpikes))
hold on
plot(T_correct, 'g')
plot(Tcorrect, 'c')