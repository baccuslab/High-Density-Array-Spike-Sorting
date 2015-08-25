pd = pdefs();
outPath    = fullfile(pd.serverData, 'HarrisBuzsaki', 'MeanShiftSorting');
figOutPath = fullfile(pd.serverData, 'HarrisBuzsaki', 'NoiseAnalysisPouzat');
if ~exist(outPath, 'file'); mkdir(outPath); end
if ~exist(figOutPath, 'file'); mkdir(figOutPath); end

anaName = 'CrossValidationThreeQuarters';

%%
for name = 1:5
    [D X Xfil I Ifil] = ana.harris.preprocess.preprocess(name);

    %%
    DS = mysort.ds.Matrix(Xfil', D.srate, 'bla');
    
    dpath = fullfile(outPath, D.name);
    [S P] = ana.sort(DS, dpath, 'r1')
    Tf = P.spikeCutting.Tf;


    %%
    vT = D.gt_template;
    nSp = size(D.gt_spikes,1);
    Eucl = sum((D.gt_spikes-repmat(vT, nSp, 1)).^2,2);
    Conv = D.gt_spikes*vT';


    %% Noise Analysis
    spikeEpochs = mysort.epoch.merge([S.spikeDetectionMerged.ts(:)-80 S.spikeDetectionMerged.ts(:)+80]);
    noiseEpochs = mysort.epoch.flip(spikeEpochs, size(DS,1));
    noiseEpochs = mysort.epoch.removeShort(noiseEpochs, P.spikeCutting.Tf);
    noiseSnippetsC = {};
    noiseEpochs_ = noiseEpochs;
    count = 0;
    while ~isempty(noiseEpochs_) && count < 50
        count = count +1;
        noiseSnippetsC{count,1} = DS.getWaveform(noiseEpochs_(:,1), 0, P.spikeCutting.Tf);
        noiseEpochs_(:,1) = noiseEpochs_(:,1)+P.spikeCutting.Tf;
        noiseEpochs_ = mysort.epoch.removeShort(noiseEpochs_, P.spikeCutting.Tf);
    end
    %%
    noiseSnippets = cell2mat(noiseSnippetsC);
    C1 = cov(noiseSnippets);
    C2 = mysort.noise.ccol2Cte(S.noise.CestS.CCol, P.spikeCutting.Tf);
        
    nNoise = size(noiseSnippets,1);
    rp = randperm(nNoise);
%     noiseSnippetsCov = noiseSnippets(floor(nNoise/2):end,:);
%     noiseSnippets = noiseSnippets(1:floor(nNoise/2)-1,:);
    noiseSnippetsCov = noiseSnippets(rp(1:floor(3*end/4)),:);
    nNoiseCov = size(noiseSnippetsCov,1);
    noiseSnippets = noiseSnippets(rp(ceil(3*end/4):end),:);
    nNoise = size(noiseSnippets,1);
    
    C = cov(noiseSnippetsCov);
    
    WN = noiseSnippets/chol(C);
    CW=cov(WN);
    MAHA = sum(WN.^2,2);
    d = size(WN,2);
    nSamples = size(WN,1);
    edges = linspace(min(MAHA)-1, max(MAHA), 50);
    binwidth = (edges(2)-edges(1));
    centers = edges(1:end-1)+binwidth/2;
    n = histc(MAHA, edges);

    %% Higher Order Moment
    nTupel = 1000;
    idx = round((d-1)*rand(nTupel,3))+1;
    HOM = zeros(nTupel,1);

    for i=1:nTupel
        HOM(i) = mean(prod(WN(:, idx(i,:)),2));    
    end
    edgeshom = linspace(min(HOM), max(HOM), 50);
    binwidthhom = (edgeshom(2)-edgeshom(1));
    centershom = edgeshom(1:end-1)+binwidthhom/2;
    nhom = histc(HOM, edgeshom);
    sdhom = 1/sqrt(nSamples);

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mysort.plot.figure([1200 800]);
    subplot(2, 4, 1:2)
    plot(centers, cumsum(n(1:end-1))/sum(n), '-k', 'linewidth', 3);
    hold on
    plot(centers, chi2cdf(centers, d-1), '-r', 'linewidth', 2)
    xlabel('Chi2')
    ylabel('Frequency')
    title(sprintf('N Samples for Cov: %d', nNoiseCov));

    subplot(2, 4, 3:4)
    bar(centers, n(1:end-1)/(binwidth*sum(n)));
    hold on
    plot(centers, chi2pdf(centers, d-1), 'g', 'linewidth', 3)
    xlabel('Chi2')
    ylabel('Density')
    title(sprintf('N Noise Samples : %d', nNoise));

    subplot(2, 4, 5:7)
    plot(HOM, 'o', 'color', [.8 .5 .5])
    hold on
    plot([1 nTupel],  3*[sdhom sdhom], ':k', 'color', [.2 .2 .2])
    plot([1 nTupel], -3*[sdhom sdhom], ':k', 'color', [.2 .2 .2])
    xlabel('Triplet Index')
    ylabel('Third Moment')

    subplot(2, 4, 8)
    barh(centershom, nhom(1:end-1)/(sum(nhom)*binwidthhom));
    hold on
    plot(normpdf(centershom, 0, sdhom), centershom, 'g', 'linewidth', 2);
    box off
    mysort.plot.figureTitle(sprintf('Harris Data set %d - %s', name, D.name));
    mysort.plot.figureDate('ana.botmpaper.harris_noise_analysis');
    mysort.plot.savefig(gcf, fullfile(figOutPath, ['noise_stat_' D.name '_' anaName]), 'fig', 0);
    % xlabel('Chi2')
    % ylabel('Density')
end
return
%%
figure
subplot(2,2,1)
imagesc(C1)
colorbar
subplot(2,2,2)
imagesc(C)
colorbar
subplot(2,2,3)
imagesc(C-C1)
colorbar

%%
figure; subplot(2,1,1)
hist(Eucl, 100)
subplot(2,1,2)
hist(Conv, 100)

%%
figure;
plot(noiseSnippets(1:10000,:)')

figure;
imagesc(C)

%%
figure;
plot(WN(1:10000,:)')
figure;
imagesc(CW)

%%
% mysort.plot.waveforms(S.clusteringMatched.templates)
mysort.plot.waveforms(S.clusteringMatched.spikeCutAligned,...
    'IDs', S.clusteringMatched.ids, 'stacked', 0)

%%
sgdf = [S.clusteringMatched.ids(:) S.clusteringMatched.ts(:)];
tgdf = [ones(length(D.gt_spiketrain),1) D.gt_spiketrain];
R = mysort.spiketrain.alignGDFs(sgdf, tgdf, 10, 10, 10);
mysort.plot.printEvaluationTable(R)

%%
for i=1:size(R.ALI,2)
    if isempty(R.ALI{i})
        continue
    end
    figure
    plot(D.gt_spikes(R.ALI{i}(:,1),:)')
end

%%
defs = mysort.util.defs();
doNotLabel = [defs.tp defs.tpo defs.fp defs.fn defs.fno];
mysort.plot.alignment(R, 'doNotLabel', doNotLabel)