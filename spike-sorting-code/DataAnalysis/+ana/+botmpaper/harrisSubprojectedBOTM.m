def = mysort.util.defs();
names = [1:4 6];
for nameidx = 4%:length(names)
    name = names(nameidx)
     fprintf('Starting with %d\n', name);
   close all
   pic = 0;
    if 1
    %     name = 2;
        [D X Xfil I Ifil] = ana.harris.preprocess.preprocess(name);

        hpf = 10; lpf= min(150, D.srate/2); forder = 4; srate = D.srate;
        hd = mysort.mea.filter_design(hpf, lpf, srate, forder);
        If                      = filtfilthd(hd, I')';     
        if 0 
            figure;
            plot(I(1:1000000))
            hold on
            plot(If(1:1000000), 'r');
        end
        if name == 1
            thr = 2200;
        elseif name == 2
            thr = 750;
        elseif name == 3
            thr = 2100;
        elseif name == 4
            thr = 2000;    
        elseif name == 6
            thr = 50;
        else
            error('threshold not set')
        end
         [pks gtSpikeTimes] = findpeaks(If+.0001*randn(size(If)), 'minpeakheight', thr, 'minpeakdistance', 10);
        if 0
            figure;
            ah = subplot(2,1,1);
            plot(I(1:2:end));
            ah(2) = subplot(2,1,2);
            plot(If(1:2:end));
            hold on
            plot(gtSpikeTimes/2, pks, 'rx', 'markersize', 20, 'linewidth', 3);
            linkaxes(ah, 'x');
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
        P.clustering.meanShiftBandWidth = sqrt(1.3*P.featureExtraction.nDims);
        P.clustering.maxSpikes = P.spikeAlignment.maxSpikes;
        
        P.botm.cutLeft = P.spikeCutting.cutLeft;
        P.botm.run = 0;
        P.botm.Tf = P.spikeCutting.Tf;
        P.botm.prior = .0001;        
        
        [S P] = ana.sort(DS, dpath, 'r2', P);
            
        Tf = P.spikeCutting.Tf;
    end
    
    %% Check feature extraction
    if 0 
        mysort.plot.figure();
        plot(S.spikePrewhitened.wfs')
        mysort.plot.spikes(S.spikePrewhitened.wfs, 'IDs', S.clustering.ids, 'nC', 4);
        mysort.plot.clusterPCA(S.spikeFeatures.X, S.clustering.ids)
%         mysort.plot.clusters(S.spikeFeatures.X, S.clustering.ids)
        mysort.plot.clustering(S.spikeFeatures.X, S.clustering.ids)
        
        
    end
    
    %%
    % Ground truth gdf:
    tgdf = [ones(length(gtSpikeTimes),1) gtSpikeTimes(:)-20];
    trueintraspikes = mysort.epoch.extractWaveform(I, [tgdf(:,2)-30 tgdf(:,2)+150]);
    trueextraspikes = DS.getWaveform(tgdf(:,2), 10, 65);
    if 0
        figure;
        subplot(2,1,1)
        plot(trueintraspikes')    
        subplot(2,1,2)
        plot(trueextraspikes')    
    end
    
    if 0 
        detTs = S.spikeDetectionMerged.ts;
        clugdf = [S.clustering.ids detTs(S.spikeAligned.alignIdx)];
        DSgt = mysort.ds.Matrix(Xfil', D.srate, 'bla');
        SC = ana.botmpaper.openSliderWithSorting(DS, clugdf, DSgt, tgdf);

        [tgdf(1:100,:) clugdf(1:100,:)]
        Rcluinv = mysort.spiketrain.alignGDFs(clugdf, tgdf, 10, 10, 10);
        rstrinv = mysort.plot.printEvaluationTable(Rcluinv);   disp(rstrinv)     
        mysort.util.logToFile(fullfile(dpath, 'rstr_msclu.txt'), ['\n' rstrinv]);
        
        Rclu = mysort.spiketrain.alignGDFs(tgdf, clugdf, 100, 10, 10);
        rstr = mysort.plot.printEvaluationTable(Rclu); disp(rstr);
        mysort.util.logToFile(fullfile(dpath, 'rstr_msclu.txt'), ['\n' rstr]);
%         mysort.util.logToFile(fullfile(dpath, 'rstr.txt'), ['\n' rstr]);        
        targetUnitID = Rclu.k2f(1)-1;
        eT = def.tp*ones(size(clugdf,1),1);
        for i=0:length(Rclu.spikeLabel2)-1
            if i==targetUnitID
                eT(clugdf(:,1)==targetUnitID) = Rclu.spikeLabel2{targetUnitID+1};
            else
                idx = Rclu.spikeLabel2{i+1}~=def.fp;
                etIdx = find(clugdf(:,1)==i);
                eT(etIdx(idx)) = Rclu.spikeLabel2{i+1}(idx);
            end
        end
        mysort.plot.clustering(S.spikeFeatures.X, S.clustering.ids, [], [], 'errorTypes', eT)
%         mysort.plot.savefig(gcf, 'ClusteringWithErrors');
        u3idx = find(clugdf(:,1)==targetUnitID);
        u3Et = Rclu.spikeLabel2{targetUnitID+1};
        mysort.plot.spikes(S.spikePrewhitened.wfs(u3idx,:), 'nC', 4, 'IDs', u3Et, 'stacked', 1);
        mysort.plot.clusterProjection(S.spikePrewhitened.wfs(u3idx,:), u3Et, [], [], 'plotTitleDiff', 1, 'plotTitleRsq', 1) 
        mysort.plot.clusterProjection(S.spikeFeatures.X(u3idx,:), u3Et, [], [], 'plotTitleDiff', 1, 'plotTitleRsq', 1) 
%         mysort.plot.savefig(gcf, 'extracellularErrors', 'fig', 0);
%         mysort.plot.spikes(S.spikePrewhitened.wfs(u3idx(u3Et==1),:), 'nC', 4, 'IDs', targetUnitID);
%         mysort.plot.spikes(S.spikePrewhitened.wfs(u3idx(u3Et==3),:), 'nC', 4, 'IDs', targetUnitID);
        wrongGtTs = clugdf(u3idx(u3Et==def.fp),2);
        intraspikes = mysort.epoch.extractWaveform(I, [wrongGtTs-100 wrongGtTs+250]);
        figure;
        plot(intraspikes')
%         mysort.plot.savefig(gcf, 'intracellularErrors', 'fig', 0);
%         figure;
%         plot(I)
%         hold on;
%         plot(wrongGtTs, I(wrongGtTs), 'rx', 'markersize', 20, 'linewidth', 4);
        
    end
    
    % Get Sorted gdf
    sgdf_ = [S.clusteringMatched.ids(:) S.clusteringMatched.ts(:)];
    templateIDs_ = unique(S.clusteringMatched.ids);
    T_ = S.clusteringMatched.templates;
    % Remove Noise unit
    if templateIDs_(1) == 0
        T = T_(2:end,:);
        templateIDs = templateIDs_(2:end);
        sgdffull = sgdf_(sgdf_(:,1)>0,:);
        cutSpikesForSorting = S.spikeCut.wfs(sgdf_(:,1)>0,:);
    else
        T = T_;
        templateIDs = templateIDs_;
        sgdffull = sgdf_;
        cutSpikesForSorting = S.spikeCut.wfs;
    end
    
    sortedUnits = unique(sgdffull(:,1));
    nSortedUnits = length(sortedUnits);
    n = histc(sgdffull(:,1), [sortedUnits-.5; sortedUnits(end)+.5]);
    n = n(1:end-1);
    unitsWithTooManySpikes = sortedUnits(find(n>4000));
    mysort.plot.figure('w', 1200, 'h', 1000);
    bar(sortedUnits, n)
    sgdf = sgdffull; %(~ismember(sgdffull(:,1), unitsWithTooManySpikes), :);
    nSp = size(D.gt_spikes,1);
    round([tgdf(1:10,2) sgdf(1:10,2)])
    Rinv = mysort.spiketrain.alignGDFs(sgdf, tgdf, 10, 10, 10);
    rstrinv = mysort.plot.printEvaluationTable(Rinv); disp(rstrinv)
    mysort.util.logToFile(fullfile(dpath, 'rstr.txt'), ['\n' rstrinv]);
    
    R = mysort.spiketrain.alignGDFs(tgdf, sgdf, 10, 10, 10);
    rstr = mysort.plot.printEvaluationTable(R); disp(rstr)
    mysort.util.logToFile(fullfile(dpath, 'rstr.txt'), ['\n' rstr]);
    
    if 0
        %% MAKE PLOT FOR MATCHING
        targetUnitID = R.k2f(1)-1;
        eT = def.tp*ones(size(sgdf_,1),1);
        for i=0:length(R.spikeLabel2)-1
            if isempty(R.spikeLabel2{i+1})
                continue
            end
            if i==targetUnitID
                eT(sgdf_(:,1)==targetUnitID+1) = R.spikeLabel2{targetUnitID+1};
            else
                idx = R.spikeLabel2{i+1}~=def.fp;
                etIdx = find(sgdf_(:,1)==i+1);
                eT(etIdx(idx)) = R.spikeLabel2{i+1}(idx);
            end
        end
        mysort.plot.clustering(S.spikeFeatures.X, S.clusteringMatched.ids(:), [], [], 'errorTypes', eT)
%         mysort.plot.savefig(gcf, 'ClusteringWithErrors');
        u3idx = find(S.clusteringMatched.ids(:)==targetUnitID);
        u3Et = R.spikeLabel2{targetUnitID+1};
        mysort.plot.spikes(S.spikePrewhitened.wfs(u3idx,:), 'nC', 4, 'IDs', u3Et, 'stacked', 1);
        mysort.plot.spikes(S.spikeCut.wfs(u3idx,:), 'nC', 4, 'IDs', u3Et, 'stacked', 1);
        mysort.plot.clusterProjection(S.spikePrewhitened.wfs(u3idx,:), u3Et, [], [], 'plotTitleDiff', 1, 'plotTitleRsq', 1) 
        mysort.plot.clusterProjection(S.spikeFeatures.X(u3idx,:), u3Et, [], [], 'plotTitleDiff', 1, 'plotTitleRsq', 1) 
%         mysort.plot.savefig(gcf, 'extracellularErrors', 'fig', 0);
%         mysort.plot.spikes(S.spikePrewhitened.wfs(u3idx(u3Et==1),:), 'nC', 4, 'IDs', targetUnitID);
%         mysort.plot.spikes(S.spikePrewhitened.wfs(u3idx(u3Et==3),:), 'nC', 4, 'IDs', targetUnitID);
        wrongGtTs = clugdf(u3idx(u3Et==def.fp),2);
        intraspikes = mysort.epoch.extractWaveform(I, [wrongGtTs-100 wrongGtTs+250]);
        figure;
        plot(intraspikes')        
    end
    
     if 0
        %% MAKE PLOT FOR FULL ONLINE BOTM
       
        sgdfBotm = S.botm.gdf;
        [tgdf(1:100,:) sgdfBotm(1:100,:)]
        sgdfBotm(:,2) = sgdfBotm(:,2)+20;
        Rbotminv = mysort.spiketrain.alignGDFs(sgdfBotm, tgdf, 10, 10, 10);
        rstrinv = mysort.plot.printEvaluationTable(Rbotminv); disp(rstrinv)
        mysort.util.logToFile(fullfile(dpath, 'botmrstr.txt'), ['\n' rstrinv]);

        Rbotm = mysort.spiketrain.alignGDFs(tgdf, sgdfBotm, 10, 10, 10);
        rstr = mysort.plot.printEvaluationTable(Rbotm); disp(rstr)
        mysort.util.logToFile(fullfile(dpath, 'botmrstr.txt'), ['\n' rstr]);        
        
        targetUnitID = Rbotm.k2f(1)-1;
        targetIdx = sgdfBotm(:,1)==targetUnitID+1;
        
        botmintraspikes = mysort.epoch.extractWaveform(I, [sgdfBotm(targetIdx,2)-50 sgdfBotm(targetIdx,2)+100]);
        botmextraspikes = DS.getWaveform(sgdfBotm(targetIdx,2), 15, 55);
        figure;
        subplot(2,1,1)
        plot(botmintraspikes')    
        subplot(2,1,2)
        plot(botmextraspikes')    
         
        
        eTBotm = def.tp*ones(size(sgdfBotm,1),1);
        for i=0:length(Rbotm.spikeLabel2)-1
            if isempty(Rbotm.spikeLabel2{i+1})
                continue
            end
            if i==targetUnitID
                eTBotm(sgdfBotm(:,1)==targetUnitID+1) = Rbotm.spikeLabel2{targetUnitID+1};
            else
                idx = Rbotm.spikeLabel2{i+1}~=def.fp;
                etIdx = find(sgdfBotm(:,1)==i+1);
                eTBotm(etIdx(idx)) = Rbotm.spikeLabel2{i+1}(idx);
            end
        end
        figure;
        subplot(2,1,1)
        plot(botmintraspikes(Rbotm.spikeLabel2{targetUnitID+1}==def.tp,:)', 'color', [.5 .5 .5])    
        hold on
        plot(botmintraspikes(Rbotm.spikeLabel2{targetUnitID+1}==def.fp,:)', 'color', [1 0 0])    
        subplot(2,1,2)
        plot(botmextraspikes(Rbotm.spikeLabel2{targetUnitID+1}==def.tp,:)', 'color', [.5 .5 .5])    
        hold on        
        plot(botmextraspikes(Rbotm.spikeLabel2{targetUnitID+1}==def.fp,:)', 'color', [1 0 0])   
        
       
    end
    
    spikeEpochs = mysort.epoch.merge([S.spikeDetectionMerged.ts(:)-120 S.spikeDetectionMerged.ts(:)+120]);
    noiseEpochs = mysort.epoch.flip(spikeEpochs, size(DS,1));
    noiseEpochs = mysort.epoch.removeShort(noiseEpochs, P.spikeCutting.Tf);
    noiseSnippets = DS.getWaveform(noiseEpochs(:,1), 0, P.spikeCutting.Tf);
    nNoise = size(noiseSnippets,1);

    XCC = mysort.noise.XCorrContainer(DS, P.spikeCutting.Tf-1, 'noiseEpochs', noiseEpochs);
    Cest = XCC.getCte4Channels(1:size(DS,2));
    dCest = diag(diag(Cest));
    



    
    % TODO: remove all correct spikes from matched cluster. replace old cluster
    % with "left over" spikes. add new cluster with correct spikes
    matchedIdx = find(R.f2k==1);
    sortedWaveformsAssignedIdx = find(sgdffull(:,1)==matchedIdx);
    matchedSpikeTrain = sgdffull(sortedWaveformsAssignedIdx,2);
    sortedWaveformsAssigned = cutSpikesForSorting(sortedWaveformsAssignedIdx,:);
    sortedWaveformsCorrectlyAssignedIdx = sortedWaveformsAssignedIdx(R.ALI{matchedIdx}(:,2));
    sortedWaveformsCorrectlyAssigned = cutSpikesForSorting(sortedWaveformsCorrectlyAssignedIdx,:);
    sortedWaveformsInCorrectlyAssignedIdx = sortedWaveformsAssignedIdx(setdiff(1:length(sortedWaveformsAssignedIdx), R.ALI{matchedIdx}(:,2)));
    sortedWaveformsInCorrectlyAssigned = cutSpikesForSorting(sortedWaveformsInCorrectlyAssignedIdx,:);
   
    correctId = find(R.f2k==1);
    correctTemplateIdx = find(templateIDs==correctId);
    correctSpikes = DS.getWaveform(D.gt_spiketrain, P.spikeCutting.cutLeft+7, P.spikeCutting.Tf);
    Tcorrect = sum(mysort.wf.v2m(mean(correctSpikes), size(DS,2)));
    [mi miidx]  = min(Tcorrect);
    [mi miidx2] = min(sum(mysort.wf.v2m(T(correctTemplateIdx,:), size(DS,2))));
    offset = miidx-miidx2;
    correctSpikes = DS.getWaveform(D.gt_spiketrain, P.spikeCutting.cutLeft+7-offset, P.spikeCutting.Tf);
    Tcorrect = mean(correctSpikes, 1);
    
    % remove spikes close to correct spikes from all spikes.
    O = mysort.spiketrain.checkForOverlaps({D.gt_spiketrain, sgdffull(:,2)}, 25);
    otherSpikes = cutSpikesForSorting(~O{2},:);
    fprintf('%d overlapping spikes removed (from %d with %d correct)\n', ...
        size(cutSpikesForSorting,1)-size(otherSpikes,1), ...
        size(cutSpikesForSorting,1), size(correctSpikes,1));
    
    

    %%
    if 0
        mysort.plot.figure('w', 1200, 'h', 1000);
        ah = subplot(3,1,1);
        plot(sortedWaveformsAssigned', 'color', [.6 .6 .6]);
        hold on
        plot(T', 'k', 'linewidth', 2)    
        plot(T(correctTemplateIdx,:), 'g', 'linewidth', 2)
        title('All assigned waveforms')
        ah(2) = subplot(3,1,2);
        plot(sortedWaveformsCorrectlyAssigned', 'color', [.6 .6 .6]);
        hold on
        plot(T', 'k', 'linewidth', 2)    
        plot(T(correctTemplateIdx,:), 'g', 'linewidth', 2)
        title('Correctly assigned waveforms')
        ah(3) = subplot(3,1,3);
        plot(sortedWaveformsInCorrectlyAssigned', 'color', [.6 .6 .6]);
        hold on
        plot(T', 'k', 'linewidth', 2)    
        plot(T(correctTemplateIdx,:), 'g', 'linewidth', 2)
        title('Incorrectly assigned waveforms')
        linkaxes(ah, 'xy');
        axis tight
        pic=pic+1; mysort.plot.savefig(gcf, fullfile(dpath, [sprintf('%03d', pic) 'assigned_and_corrAssigned_waveforms']), 'fig', 0);

        mysort.plot.figure('w', 1200, 'h', 1000);

        plot(cutSpikesForSorting', 'color', [.6 .6 .6]);
        hold on
        plot(T', 'k', 'linewidth', 2)    
        plot(T(correctTemplateIdx,:), 'g', 'linewidth', 2)
        pic=pic+1; mysort.plot.savefig(gcf, fullfile(dpath, [sprintf('%03d', pic) 'corr_template_over_all_templates']), 'fig', 0);

        mysort.plot.figure('w', 1200, 'h', 1000);

        plot(correctSpikes', 'k')
        hold on
        plot(T(correctTemplateIdx,:), 'g', 'linewidth', 2)
        plot(Tcorrect, 'b', 'linewidth', 2)
        pic=pic+1; mysort.plot.savefig(gcf, fullfile(dpath, [sprintf('%03d', pic) 'corr_template_over_corr_spikes']));

        mysort.plot.figure('w', 1200, 'h', 1000);

        plot(T');
        plot(T(correctTemplateIdx,:), 'g', 'linewidth', 2)
        hold on
        plot(D.gt_template, 'k', 'linewidth', 2)
        pic=pic+1; mysort.plot.savefig(gcf, fullfile(dpath, [sprintf('%03d_', pic) 'all_templates']));
    end

    %%
%     Eucl  = @(X, t)    sum(  (X-repmat(t, size(X,1), 1)).^2                                ,2);
%     Maha  = @(X, t, C) sum( ((X-repmat(t, size(X,1), 1))/C) .* (X-repmat(t, size(X,1), 1)) ,2);
%     Conv  = @(X, t)    X*t';
%     Match = @(X, t, C) X*(t/C)';
%     Botm  = @(X, t, C) X*(t/C)' - .5*t*(t/C)';
    Eucl  = @ana.botmpaper.Eucl;
    Maha  = @ana.botmpaper.Maha;
    Conv  = @ana.botmpaper.Conv;
    Match = @ana.botmpaper.Match;
    Botm  = @ana.botmpaper.Botm;

    C = cov(noiseSnippets);
    dC = diag(diag(C));
    targetC = dC;
    targetCest = dCest;
    %% see Pope 2008, "Shrinkage estimation of the power spectrum covariance matrix", eq.10
    n = size(noiseSnippets,1);
    meanNoise = mean(noiseSnippets,1);
    meanFreeNoise = noiseSnippets - repmat(meanNoise,n,1);
    varC  = zeros(size(C));
    covTC = varC;
%     varCest  = zeros(size(Cest));
%     covTCest = varCest;
    for k=1:n
        Wijk = meanFreeNoise(k,:)'*meanFreeNoise(k,:);
        Tijk = diag(diag(Wijk));
        varC = varC + (Wijk - C).^2;
        covTC = covTC + (Wijk - C).*(Tijk - targetC);
    end
    varC = varC * n/(n-1)^3;
    covTC = covTC * n/(n-1)^3;

    %%
    lambdaopt = sum(varC(:) - covTC(:))   /  sum(sum((C-targetC)).^2);
    lambda = .5;
    CDL = ana.botmpaper.diagonalLoading(C, 10000);
    C_ = (1-lambda)*C + lambda*targetC;
    C_opt = (1-lambdaopt)*C + lambdaopt*targetC;
    Cest_ = (1-lambda)*Cest + lambda*targetCest;
    CestDL = ana.botmpaper.diagonalLoading(Cest, 10000);

    %%   SORT WITH BOTM ON FULL DATA
    if 0
        DS1 = mysort.ds.Matrix(Xfil', D.srate, 'bla');
        MSSortingContainer = mysort.spiketrain.SpikeSortingContainer('meanShift', sgdf, 'wfDataSource', DS);
        DS1.addSpikeSorting(MSSortingContainer);    

        DS2 = mysort.ds.Matrix(Xfil', D.srate, 'bla');
        TrueSortingContainer = mysort.spiketrain.SpikeSortingContainer('gt', tgdf, 'wfDataSource', DS);
        DS2.addSpikeSorting(TrueSortingContainer);

        Covest = struct();
        Covest.CCol = mysort.noise.Cte2Ccol(C_, size(DS,2));
        T_ = mysort.wf.vSubsel(T,size(DS,2),2:(size(T,2)/size(DS,2))-1);
    
        botm = mysort.sorters.BOTMOvp(Covest, mysort.wf.v2t(T_, size(DS,2)), 'adaptOnInit', 1);
        botm.DH = DS;
        botm.sortMatrix(DS(1:1000,:)');
        mysort.plot.SliderDataAxes({DS1, DS2, botm}, 'channelSpacers', [400, 400, 0])
        % samples : 2.42652e+006

        botmgdf1 = botm.sort(DS);
    %     botmgdf2 = botm.sort(DS);
    %     save('b2_minusShift_gdf', 'botmgdf2')
    %     load('b1_plusShift_gdf.mat', 'botmgdf1'); 
    %     load('b2_minusShift_gdf.mat', 'botmgdf2'); 

        botmgdf1(:,2) = botmgdf1(:,2)-ceil((size(T,2)/size(DS,2))/2)-1;
    %     botmgdf2(:,2) = botmgdf2(:,2)-ceil((size(T,2)/size(DS,2))/2)-1;

        DS3 = mysort.ds.Matrix(Xfil', D.srate, 'bla');
        BOTMSortingContainer1 = mysort.spiketrain.SpikeSortingContainer('gt', botmgdf1, 'wfDataSource', DS);
        DS3.addSpikeSorting(BOTMSortingContainer1);

    %     DS4 = mysort.ds.Matrix(Xfil', D.srate, 'bla');
    %     BOTMSortingContainer2 = mysort.spiketrain.SpikeSortingContainer('gt', botmgdf2, 'wfDataSource', DS);
    %     DS4.addSpikeSorting(BOTMSortingContainer2);

        mysort.plot.SliderDataAxes({DS1, DS2, DS3, botm}, 'channelSpacers', [400, 400, 400, 0])
    %     mysort.plot.SliderDataAxes({DS1, DS2, DS3, DS4, botm}, 'channelSpacers', [400, 400, 400, 400, 0])
        botmgdf1_ = botmgdf1(botmgdf1(:,1)>1,:);
        tgdf_ = tgdf;
        tgdf_(:,2) = tgdf_(:,2)-23;
        [tgdf_(1:10,:) botmgdf1_(1:10,:)]

        [botmgdf1__ wasResolved] = mysort.spiketrain.resolveOvpIndexInGdf(botmgdf1_, botm.DOvpC.DOvpIndex);
        Rinv = mysort.spiketrain.alignGDFs(botmgdf1__, tgdf_, 10, 50, 10);
        rstrinv = mysort.plot.printEvaluationTable(Rinv); disp(rstrinv)
        mysort.util.logToFile(fullfile(dpath, 'rstr_botm.txt'), ['\n' rstrinv]);

        R = mysort.spiketrain.alignGDFs(tgdf_, botmgdf1__, 10, 50, 10);
        rstr = mysort.plot.printEvaluationTable(R); disp(rstr)
        mysort.util.logToFile(fullfile(dpath, 'rstr_botm.txt'), ['\n' rstr]);    
        spikeLabels = R.spikeLabel2{R.St2IDs==3};
        spikeLabelsWasResolved = wasResolved(botmgdf1__(:,1)==3);
        unit3wfs = DS.getWaveform(botmgdf1__(botmgdf1__(:,1)==3,2), 10, 60);
        %%
        fpIdx = ismember(spikeLabels, [def.fp def.fpo]);
        tpIdx = ismember(spikeLabels, [def.tp def.tpo]);
        if 0
            figure;
            ah = subplot(5,1,1);
            plot(unit3wfs')
            title(['all spikes (' num2str(size(unit3wfs,1)) ')'])

            ah(2) = subplot(5,1,2);
            Sp = unit3wfs(tpIdx' & ~spikeLabelsWasResolved,:)';
            plot(Sp)
            title(['tp spikes not resolved (' num2str(size(Sp,2)) ')'])

            ah(3) = subplot(5,1,3);
            Sp = unit3wfs(tpIdx' & spikeLabelsWasResolved,:)';
            plot(Sp)
            title(['tp spikes resolved (' num2str(size(Sp,2)) ')'])   

            ah(4) = subplot(5,1,4);
            Sp = unit3wfs(fpIdx' & ~spikeLabelsWasResolved,:)';
            plot(Sp)
            title(['fp spikes not resolved (' num2str(size(Sp,2)) ')'])    

            ah(5) = subplot(5,1,5);
            Sp = unit3wfs(fpIdx' & spikeLabelsWasResolved,:)';
            plot(Sp)
            title(['fp spikes resolved (' num2str(size(Sp,2)) ')'])    
            linkaxes(ah, 'xy');   
        end
    end
    %%
    if 0
        mysort.plot.figure('w', 1200, 'h', 1000);

        subplot(2,2,1)
        imagesc(C)
        colorbar
        title('C')

        subplot(2,2,2)
        imagesc(C_opt)
        colorbar
        title('Cloaded')

        subplot(2,2,3)
        imagesc(varC)
        colorbar
        title('var(Cij)')

        subplot(2,2,4)
        imagesc(covTC)
        colorbar
        title('cov(Tij,Cij)')
        pic=pic+1; mysort.plot.savefig(gcf, fullfile(dpath, [sprintf('%03d_', pic) 'covs']));
    end
    %%
    H = cov(cutSpikesForSorting);

    [VH DH] = eig(H);
    [VC DC] = eig(C);
    [VC_ DC_] = eig(C_);
    dH = diag(DH);
    largeEV = dH > .001*max(dH);
    Proj = VH(:, largeEV);
    fprintf('Chosen number of dimensions: %d\n', size(Proj,2));

    CP = cov(noiseSnippets*Proj);
    CPest = Proj'*Cest*Proj;
    HP = cov([otherSpikes; correctSpikes]);
    TP = T*Proj;

    % This weird way of passing the arguments became necessary to estimate
    % the computation times correctly. Some preprocessing can be done
    % before the template matching is started, some must be done later
    DATASETS = {otherSpikes, T, C, 'other, C', 0, []
                correctSpikes,  T, C, 'cor, C', 0, []
                noiseSnippets,  T, C, 'N, C', 0, []
                otherSpikes, T, dC, 'other, dC', 0, []
                correctSpikes,  T, dC, 'cor, dC', 0, []
                noiseSnippets,  T, dC, 'N, dC', 0, []
                otherSpikes, T, C_opt, 'other, C_opt', 0, []
                correctSpikes,  T, C_opt, 'cor, C_opt', 0, []
                noiseSnippets,  T, C_opt, 'N, C_opt', 0, []
                otherSpikes, T, C_, 'other, C_', 0, []
                correctSpikes,  T, C_, 'cor, C_', 0, []
                noiseSnippets,  T, C_, 'N, C_' , 0  , []             
                otherSpikes, TP, CP, 'other, CP', 0, Proj
                correctSpikes, TP, CP, 'cor, CP', 0, Proj
                noiseSnippets, TP, CP, 'N, CP' , 0 , Proj   
                
                otherSpikes, T, Cest, 'other, Cest', 0, []
                correctSpikes,  T, Cest, 'cor, Cest', 0, []
                noiseSnippets,  T, Cest, 'N, Cest', 0, []
                otherSpikes, T, dCest, 'other, dC', 0, []
                correctSpikes,  T, dCest, 'cor, dC', 0, []
                noiseSnippets,  T, dCest, 'N, dC' , 0 , [] 
                otherSpikes, T, Cest_, 'other, Cest_', 0, []
                correctSpikes,  T, Cest_, 'cor, Cest_', 0, []
                noiseSnippets,  T, Cest_, 'N, Cest_', 0, []
                otherSpikes, TP, CPest, 'other, CPest', 0, Proj
                correctSpikes, TP, CPest, 'cor, CPest', 0, Proj
                noiseSnippets, TP, CPest, 'N, CPest', 0, Proj
                otherSpikes, T, CDL, 'other, CDL', 0, []
                correctSpikes,  T, CDL, 'cor, CDL', 0, []
                noiseSnippets,  T, CDL, 'N, CDL', 0, []
                otherSpikes, T, CestDL, 'other, CDL', 0, []
                correctSpikes,  T, CestDL, 'corr, CestDL', 0, []
                noiseSnippets,  T, CestDL, 'N, CestDL', 0, []};
            
    datasetNames = DATASETS(:,4);
        
    methods = { Eucl, 'noC', 'minDetection', 'Eucl'
                Maha, 'C',   'minDetection', 'Maha'
                Conv, 'noC', 'maxDetection', 'Conv'
                Match, 'C',  'maxDetection', 'Match'
                Botm, 'C',   'maxDetection', 'Botm'
                Botm, 'C',   'maxDetection', 'BotmOpt'};      
            % for BotmOpt the detectin threshold will be set to 0 (close to
            % the theoretical optimum) for 'Botm' is is calculated by
            % fminsearch (see below stcmp ... 'Opt' )

    %%
    mysort.plot.figure('w', 1200, 'h', 1000);
    
    plot(diag(DH))
    hold on
    plot(diag(DC), 'r')
    plot(diag(DC_), 'r:')
    legend('H', 'C', 'C.')
    title('Eigenvalues of covariance matrices');
    pic=pic+1; mysort.plot.savefig(gcf, fullfile(dpath, [sprintf('%03d_', pic) 'eigenvalues_of_covs']));

    %% Compute Classifications
    nD = size(DATASETS,1);
    nMethods = size(methods,1);
    RES = [];
    nBins = 400;
    for d=1:nD
        fprintf('Dataset: %d\n', d);
        nSp = size(DATASETS{d,1},1);
        lT = DATASETS{d,2};
        nT  = size(lT,1);
        RES(d).TM = zeros(nSp, nMethods, nT);
        RES(d).counts = zeros(nBins, nMethods, nT);
        for m=1:nMethods
            fprintf('Method: %d\n', m);
            for i=1:nT
                fun = methods{m,1};
                t1 = tic;
                if strcmp(methods{m,2}, 'C')
                    RES(d).TM(:,m,i) = fun(DATASETS{d,1}, lT(i,:), DATASETS{d,3}, DATASETS{d,6});
                    RES(d).cputime(m,i) = toc(t1)+DATASETS{d,5};
%                     RES(d).templateEnergies(i) = sqrt((lT(i,:)/DATASETS{d,3})*lT(correctTemplateIdx,:)');
                else
                    RES(d).TM(:,m,i) = fun(DATASETS{d,1}, lT(i,:), DATASETS{d,6});
                    RES(d).cputime(m,i) = toc(t1)+DATASETS{d,5};
%                     RES(d).templateEnergies(i) = norm(lT(i,:));
                end
            end
            fprintf('Estimating performance...\n');
            if strcmp(methods{m,3}, 'minDetection')
                [mx mxidx] = min(squeeze(RES(d).TM(:,m,:)),[],2);
            else
                [mx mxidx] = max(squeeze(RES(d).TM(:,m,:)),[],2);
            end            
            RES(d).classifications(:,m) = histc(mxidx, [0.5:1:(nT+0.5)]);
            RES(d).performance(m) = sum(mxidx==correctTemplateIdx)/length(mxidx);   

            fprintf('Estimating bins...\n');
            RES(d).minPerMethod(m) = min(min(RES(d).TM(:,m,:)));
            RES(d).maxPerMethod(m) = max(max(RES(d).TM(:,m,:)));
            RES(d).edges(:,m) = linspace(RES(d).minPerMethod(m), RES(d).maxPerMethod(m), nBins);
            RES(d).widthPerMethod(m) = RES(d).edges(2,m) - RES(d).edges(1,m);
            RES(d).binCenters(:,m) = RES(d).edges(:,m) + RES(d).widthPerMethod(m)*.5;
            fprintf('Counting...\n');
            RES(d).counts(:,m,:) = histc(squeeze(RES(d).TM(:,m,:)), RES(d).edges(:,m));
        end
    end
    fprintf('Done.\n');

    %% Make plot of count histograms for noise and signal
    noiseCorrectPairs = repmat([2  3], size(DATASETS,1)/3, 1) + repmat((0:(size(DATASETS,1)/3) -1)'*3, 1, 2);
    detThr = zeros(size(noiseCorrectPairs,1), nMethods);
    detPerf = detThr;
    mysort.plot.figure('w', 1200, 'h', 1000);
    p = 0;
    for pairi = 1:size(noiseCorrectPairs,1)
        correctDSIdx = noiseCorrectPairs(pairi,1);
        noiseDSIdx = noiseCorrectPairs(pairi,2);

        neuronFR = 30; %Hz
        samplesPerSecond = 10000;
        noisePriorWeight = 1; %samplesPerSecond/neuronFR;
        nCorrect = size(RES(correctDSIdx).TM,1);
        nNoise   = size(RES(noiseDSIdx).TM,1);

        for m=1:nMethods
            p = p+1;
            subplot(size(noiseCorrectPairs,1), nMethods, p);
            plot(RES(noiseDSIdx).binCenters(:,m), RES(noiseDSIdx).counts(:,m,correctTemplateIdx) , '.-r')
            hold on
            plot(RES(correctDSIdx).binCenters(:,m), RES(correctDSIdx).counts(:,m,correctTemplateIdx) , '.-g')
            miidx = find(RES(correctDSIdx).counts(:,m,correctTemplateIdx)>2,1);
            maidx = length(RES(correctDSIdx).counts(:,m,correctTemplateIdx)) - find(flipud(RES(correctDSIdx).counts(:,m,correctTemplateIdx))>2,1);
            if maidx-miidx < 20
                miidx = max(1, miidx-20);
                maidx = min(length(RES(correctDSIdx).counts(:,m,correctTemplateIdx)), maidx+20);
            end
            if m==1
                ylabel(DATASETS{noiseDSIdx, 4});
            end
            if pairi ==1
                title(methods{m, 4});
            end
            set(gca, 'xlim', sort([RES(correctDSIdx).binCenters(miidx,m) RES(correctDSIdx).binCenters(maidx,m)]));

        end
    end
    pic=pic+1; mysort.plot.savefig(gcf, fullfile(dpath, [sprintf('%03d_', pic) 'counts_noise_signal']));
    %% Use Classification responses to noise and to the correct template to compute the
    % detection performance

    detThr = zeros(size(noiseCorrectPairs,1), nMethods);
    detPerf = detThr;
    mysort.plot.figure('w', 1200, 'h', 1000);
    p = 0;

    for pairi = 1:size(noiseCorrectPairs,1)
        correctDSIdx = noiseCorrectPairs(pairi,1);
        noiseDSIdx = noiseCorrectPairs(pairi,2);

        neuronFR = 30; %Hz
        samplesPerSecond = 10000;
        noisePriorWeight = samplesPerSecond/neuronFR;
        nCorrect = size(RES(correctDSIdx).TM,1);
        nNoise   = size(RES(noiseDSIdx).TM,1);

        for m=1:nMethods
            noisefun = @(thr) sum(RES(noiseDSIdx).TM(:,m,correctTemplateIdx) < thr)/nNoise;
            signafun = @(thr) sum(RES(correctDSIdx).TM(:,m,correctTemplateIdx) < thr)/nCorrect;

            if strcmp(methods{m,3}, 'minDetection')
                E = @(thr)  (noisePriorWeight*noisefun(thr) - signafun(thr))/noisePriorWeight;
                thr0 = .8 * mean(RES(noiseDSIdx).TM(:,m,correctTemplateIdx));
            else
                E = @(thr) (-noisePriorWeight*noisefun(thr) + signafun(thr))/noisePriorWeight;        
                thr0 = .8 * mean(RES(correctDSIdx).TM(:,m,correctTemplateIdx));
            end
            if ~isempty(strfind(methods{m, 4}, 'Opt'))
                detThr(pairi,m) = 0; 
                detPerf(pairi, m) = E(0);
            else            
                [detThr(pairi,m) detPerf(pairi, m)] = fminsearch(E, thr0);
            end
            if 1
                mima = sort([-1*detThr(pairi,m) 3*detThr(pairi,m)]);
                mima = [min(-50, mima(1)) max(100, mima(2))];
                xrange = linspace(mima(1),mima(2), 1000); %0:.05:10;
                K = zeros(1, length(xrange));
                KS = K;
                KN = K;
                for i=1:length(xrange)
                    KS(i) = signafun(xrange(i));
                    KN(i) = noisefun(xrange(i));
                    K(i) = E(xrange(i));
                end
                p = p+1;
                subplot(size(noiseCorrectPairs,1), nMethods, p);
                plot(xrange, K, '.-')
                hold on
                plot(xrange, KS, '.-g')
                plot(xrange, KN, '.-r')
                plot(thr0, 0, 'bd');
                plot(detThr(pairi,m), detPerf(pairi, m), 'gd');
                if m==1
                    ylabel(DATASETS{noiseDSIdx, 4});
                end
                if pairi ==1
                    title(methods{m, 4});
                end
            end
        end
    end

    pic=pic+1; mysort.plot.savefig(gcf, fullfile(dpath, [sprintf('%03d_', pic) 'detection_performance_estimation']));

    %% Make big classification plot, all vs. all
    mysort.plot.figure('w', 1200, 'h', 1000);
    
    p = 1;
    ah = [];
    correctSpikesIdx = 2:3:size(DATASETS,1);
    for d=correctSpikesIdx
        lT = DATASETS{d,2};
        nT  = size(lT,1); 
        nS  = size(DATASETS{d,1},1);
        for m=1:nMethods
            ah(p) = subplot(length(correctSpikesIdx),nMethods,p);
            bar(1:nT, RES(d).classifications(1:nT,m))
            title(methods{m,4});
            if m == 1
                ylabel(DATASETS{d, 4})
            end
            p=p+1;
        end
    end
    set(ah, 'ylim', [0 nS+10]);
    mysort.plot.figureTitle('Correct Classification');
%     pic=pic+1; mysort.plot.savefig(gcf, fullfile(dpath, [sprintf('%03d_', pic) 'classifications']));
    %% Make big plot of all errors, all vs. all
    mysort.plot.figure('w', 1200, 'h', 1000);
    
    p = 1;
    ah = [];
    correctSpikesIdx = 2:3:size(DATASETS,1);
    for d=correctSpikesIdx
        lT = DATASETS{d,2};
        nT  = size(lT,1); 
        nS  = size(DATASETS{d,1},1);
        for m=1:nMethods
            ah(p) = subplot(length(correctSpikesIdx),nMethods,p);
            bar(1:nT, RES(d).classifications(1:nT,m))
            title(methods{m,4});
            if m == 1
                ylabel(DATASETS{d, 4})
            end
            p=p+1;
        end
    end
    set(ah, 'ylim', [0 nS+10]);
    mysort.plot.figureTitle('Correct Classification');
    pic=pic+1; mysort.plot.savefig(gcf, fullfile(dpath, [sprintf('%03d_', pic) 'all errors']));    
    %% Make big performance plot, only performance
    ah = [];
    PERF = [];
    for di=1:length(correctSpikesIdx)
        d = correctSpikesIdx(di);
        for m=1:nMethods
            PERF(di,m) = 100*RES(d).classifications(correctTemplateIdx,m)/sum(RES(d).classifications(:,m));
        end
    end
    otherSpikesIdx = 1:3:(size(DATASETS,1)-1);
    REJ = [];
    for di=1:length(otherSpikesIdx)
        d = otherSpikesIdx(di);
        for m=1:nMethods
            REJ(di,m) = 100 - 100*RES(d).classifications(correctTemplateIdx,m)/sum(RES(d).classifications(:,m));
        end
    end    
    mysort.plot.figure('w',1100, 'h', 800);
    ah = subplot(2,2,1);
    ah(2) = subplot(2,2,2);
    ah(3) = subplot(2,2,3);
    ah(4) = subplot(2,2,4);    
    set(ah, 'fontsize', 12);
    set(ah(4), 'xticklabel', DATASETS(correctSpikesIdx,4))
    
    bar(ah(1), -100*detPerf)
    ylabel(ah(1), 'Detection Perf [%]', 'fontsize', 12);
    set(ah(1), 'ylim', [60 100]);    
    
    bar(ah(2), PERF);    
    set(ah(2), 'ylim', [0 100]);
    ylabel(ah(2), 'Class Perf [%]', 'fontsize', 12);

    bar(ah(3), REJ);    
    set(ah(3), 'ylim', [0 100]);
    ylabel(ah(3), 'Reject Perf [%]', 'fontsize', 12);    
    
    tot = (-100*detPerf + PERF + REJ)/3;
    bar(ah(4), tot);    
    set(ah(4), 'ylim', [0 100]);
    ylabel(ah(4), 'Total [%]', 'fontsize', 12);  
    set(ah(4), 'ylim', [min(tot(:)) 100]);
    
    legend(ah(1), methods(:,4), 'location', 'northeastoutside')
    legend(ah(2), methods(:,4), 'location', 'northeastoutside')    
    legend(ah(3), methods(:,4), 'location', 'northeastoutside')
    
    pic=pic+1; mysort.plot.savefig(gcf, fullfile(dpath, [sprintf('%03d_', pic) 'final_results']));
    
    %% Make plot of theoretical and actual distribution of norms
    mysort.plot.figure('w', 1200, 'h', 1000);
    p = 1;
    nP = size(noiseCorrectPairs,1);
    ah = [];    
    for pairi = 1:nP
        correctDSIdx = noiseCorrectPairs(pairi,1);
        noiseDSIdx = noiseCorrectPairs(pairi,2);
        S      = DATASETS{correctDSIdx,1};
        N      = DATASETS{noiseDSIdx,1};
        if ~isempty(DATASETS{correctDSIdx,6})
            S = S*DATASETS{correctDSIdx,6};
            N = N*DATASETS{correctDSIdx,6};
        end
        lT     = DATASETS{correctDSIdx,2};
        Cpair  = DATASETS{correctDSIdx,3};
        U = chol(Cpair);
        
%         Cnoise = DATASETS{noiseDSIdx,3};
        
        dims = size(S,2);
        corrT = repmat(mean(S), size(S,1), 1);
        R = S - corrT;
        R = R/U;
        e = norm(corrT/U);
        N = N/U;
        
        correct_distr = sum(R.^2,2);
%         correct_edges = linspace(min(correct_distr), max(correct_distr), 1000);
        if pairi == nP
            correct_edges = linspace(0, 110, 100);
        else
            correct_edges = linspace(0, 400, 100);
        end
        wcorr = correct_edges(2) - correct_edges(1);
        correct_bincenters = correct_edges + wcorr/2;
        correct_distr_counts = histc(correct_distr, correct_edges);
        correct_distr_pdf    = correct_distr_counts/(sum(correct_distr_counts)*wcorr);
        theo_correct = chi2pdf((correct_edges), dims);
        
        noise_distr = sum(N.^2,2);
%         noise_edges = linspace(min(noise_distr), max(noise_distr), 1000);
        noise_edges = correct_edges;
        wnoise = noise_edges(2) - noise_edges(1);
        noise_bincenters = noise_edges + wnoise/2;        
        noise_distr_counts = histc(noise_distr, noise_edges);
        noise_distr_pdf    = noise_distr_counts/(sum(noise_distr_counts)*wnoise);
        theo_noise  = theo_correct; %ncx2pdf(noise_edges, dims, e^2);
        
        ah(p) = subplot(2,nP,p); hold on
        plot(ah(p), correct_bincenters, correct_distr_pdf, 'g', 'linewidth', 2);
        plot(ah(p), correct_bincenters, theo_correct, 'k', 'linewidth', 2);
        title(DATASETS{correctDSIdx,4});
        
        ah(p) = subplot(2,nP,p+nP); hold on
        plot(ah(p), noise_bincenters, noise_distr_pdf, 'r', 'linewidth', 2);
        plot(ah(p), noise_bincenters, theo_noise, 'k', 'linewidth', 2);
        axis tight
        p = p+1;
    end
%     pic=pic+1; mysort.plot.savefig(gcf, fullfile(dpath, [sprintf('%03d_', pic) 'distributionComparisons']));
    
% scalefactor = 1/mean(sqrt(diag(C)));
% scalefactor = 1/max((sqrt(diag(C))));
% mysort.plot.figure('w', 1000, 'h', 1000);
% p=0; nR = 5; nC=2;
% % Eucl
% p=p+1; ah(p) = subplot(nR,nC,p);
% dims = size(T,2);
% width = edges(2)-edges(1);
% scalefactor = 1/mean(sqrt(diag(C_)));
% % Eucl
% % plot(ah(p), edges, nEucl(:,correctIdx)/(sum(nEucl(:,correctIdx))*width), 'b', 'linewidth', 2);
% hold on
% plot(ah(p), edges, nEuclCor(:,correctIdx)/(sum(nEuclCor(:,correctIdx))*width), 'g', 'linewidth', 2);
% e = norm(T(correctIdx,:));
% % theo = 10e4*scalefactor^2*ncx2pdf(scalefactor^2*(edges), Tf, (scalefactor*e)^2);
% 
% plot(ah(p), edges/10, theo_correct, 'g:', 'linewidth', 2);
% plot(ah(p), edges/10, theo_noise, 'r:', 'linewidth', 2);
% plot(ah(p), edges, nEuclNoise(:,correctIdx)/(sum(nEuclNoise(:,correctIdx))*width), 'r', 'linewidth', 2);
% E = @(thr) sum(nEuclCor(:,correctIdx)<thr) 
% p=p+1; ah(p) = subplot(nR,nC,p);
% [mx mxidx] = min(EuclCor,[],2);
% n = histc(mxidx, [0:1:12]);
% classPerfEucl = sum(mxidx==correctIdx)/length(mxidx);
% bar(1:12, n(1:12))
% title('Eucl');    
    
    
    
    %%
    clear X Xfil I Ifil S DATASETS noiseSnippets XCC D meanFreeNoise SP NP mx mxidx K KS KN otherSpikes
%     save(fullfile(dpath, 'final_results.mat'))
    save(fullfile(dpath, 'final_results_small.mat'), 'correctTemplateIdx', 'RES', 'detPerf', 'methods', 'datasetNames');
    fprintf('Done with %d\n', name);
end