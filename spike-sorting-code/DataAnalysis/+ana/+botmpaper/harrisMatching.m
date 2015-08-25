if 0 
    name = 1;
[D X Xfil I Ifil] = ana.harris.preprocess.preprocess(name);

%%
DS  = mysort.ds.Matrix(Xfil', D.srate, 'Prefiltered Extracellular');
DSI = mysort.ds.Matrix(I', D.srate,    'Raw Intracellular');
% dpath = ['C:\LocalData\HarrisBuzsaki\MeanShiftSorting\' D.folder];
dpath = ana.harris.datapath();
dpath = fullfile(dpath, '..', 'MeanShiftSorting6', D.name);
[S P] = ana.sort(DS, dpath, 'r1')
Tf = P.spikeCutting.Tf;
%%
% mysort.plot.waveforms(S.clusteringMatched.templates)
mysort.plot.waveforms(S.clusteringMatched.spikeCutAligned,...
    'IDs', S.clusteringMatched.ids, 'stacked', 0)
end
%%
sgdf = [S.clusteringMatched.ids(:) S.clusteringMatched.ts(:)];
tgdf = [ones(length(D.gt_spiketrain),1) D.gt_spiketrain];
R = mysort.spiketrain.alignGDFs(sgdf, tgdf, 10, 10, 10);
mysort.plot.printEvaluationTable(R)

%%
defs = mysort.util.defs();
doNotLabel = [defs.tp defs.tpo defs.fp defs.fn defs.fno];
mysort.plot.alignment(R, 'doNotLabel', doNotLabel)

%%
for i=1:size(R.ALI,2)
    if isempty(R.ALI{i})
        continue
    end
    figure
    plot(D.gt_spikes(R.ALI{i}(:,1),:)')
end

%%
vT = D.gt_template;
nSp = size(D.gt_spikes,1);
Eucl = sum((D.gt_spikes-repmat(vT, nSp, 1)).^2,2);
Conv = D.gt_spikes*vT';

%%
figure; subplot(2,1,1)
hist(Eucl, 100)
subplot(2,1,2)
hist(Conv, 100)
%% Noise Analysis
spikeEpochs = mysort.epoch.merge([S.spikeDetectionMerged.ts(:)-80 S.spikeDetectionMerged.ts(:)+80]);
noiseEpochs = mysort.epoch.flip(spikeEpochs, size(DS,1));
noiseEpochs = mysort.epoch.removeShort(noiseEpochs, P.spikeCutting.Tf);
noiseSnippets = DS.getWaveform(noiseEpochs(:,1), 0, P.spikeCutting.Tf);
nNoise = size(noiseSnippets,1);

%%
figure;
plot(noiseSnippets(1:min(end, 10000),:)')
C=cov(noiseSnippets);
figure;
imagesc(C)

%%
T = S.clusteringMatched.templates;
correctId = R.f2k;
correctSpikes = DS.getWaveform(D.gt_spiketrain, P.spikeCutting.cutLeft+7, P.spikeCutting.Tf);

nSp = size(S.spikeCut.wfs,1);
nSpCor = size(correctSpikes,1);
nT = size(T,1);

Eucl = [];
EuclCor = [];
EuclNoise = [];

Maha = [];
MahaCor = [];
MahaNoise = [];

Conv = [];
ConvCor = [];
ConvNoise = [];

Matched = [];
MatchedCor = [];
MatchedNoise = [];

botm = [];
botmCor = [];
botmNoise = [];

nEucl = [];
nEuclCor = [];
nEuclNoise = [];

nConv = [];
nConvCor = [];
nConvNoise = [];

nBotm = [];
nBotmCor = [];
nBotmNoise = [];
C_ = .5*C + .5*diag(diag(C));
U_ = chol(C_);
for i=1:nT
    F = C_\T(i,:)';
        
    Eucl(:,i)    = sum((S.spikeCut.wfs-repmat(T(i,:), nSp, 1)).^2,2);
    EuclCor(:,i) = sum((correctSpikes-repmat(T(i,:), nSpCor, 1)).^2,2);
    EuclNoise(:,i) = sum((noiseSnippets-repmat(T(i,:), nNoise, 1)).^2,2);

    Maha(:,i)    = sum( ((S.spikeCut.wfs-repmat(T(i,:), nSp, 1))/C_).*(S.spikeCut.wfs-repmat(T(i,:), nSp, 1)) ,2);
    MahaCor(:,i) = sum( ((correctSpikes-repmat(T(i,:), nSpCor, 1))/C_).*(correctSpikes-repmat(T(i,:), nSpCor, 1)) ,2);
    MahaNoise(:,i) = sum( ((noiseSnippets-repmat(T(i,:), nNoise, 1))/C_).*(noiseSnippets-repmat(T(i,:), nNoise, 1)) ,2);

    Conv(:,i)    = S.spikeCut.wfs*T(i,:)';
    ConvCor(:,i) = correctSpikes*T(i,:)';
    ConvNoise(:,i) = noiseSnippets*T(i,:)';
    
    Matched(:,i)    = S.spikeCut.wfs*F;
    MatchedCor(:,i) = correctSpikes*F;
    MatchedNoise(:,i) = noiseSnippets*F;    
    
    botm(:,i)    = S.spikeCut.wfs*F - .5*T(i,:)*F;
    botmCor(:,i) = correctSpikes*F  - .5*T(i,:)*F;
    botmNoise(:,i) = noiseSnippets*F  - .5*T(i,:)*F;
    
    edges = linspace(-10e8, 10e8, 20000);
    botm_edges = linspace(-10e4, 10e4, 20000);
    
    nEucl(:,i) = histc(Eucl(:,i), edges);
    nEuclCor(:,i) = histc(EuclCor(:,i), edges);
    nEuclNoise(:,i) = histc(EuclNoise(:,i), edges);
    
    nMaha(:,i) = histc(Maha(:,i), botm_edges);
    nMahaCor(:,i) = histc(MahaCor(:,i), botm_edges);
    nMahaNoise(:,i) = histc(MahaNoise(:,i), botm_edges);
    
    nConv(:,i) = histc(Conv(:,i), edges);
    nConvCor(:,i) = histc(ConvCor(:,i), edges);
    nConvNoise(:,i) = histc(ConvNoise(:,i), edges);    
    
    nMatched(:,i) = histc(Matched(:,i), botm_edges);
    nMatchedCor(:,i) = histc(MatchedCor(:,i), botm_edges);
    nMatchedNoise(:,i) = histc(MatchedNoise(:,i), botm_edges);
    
    nBotm(:,i) = histc(botm(:,i), botm_edges);
    nBotmCor(:,i) = histc(botmCor(:,i), botm_edges);
    nBotmNoise(:,i) = histc(botmNoise(:,i), botm_edges);
end

%%
DS = mysort.ds.Matrix(Xfil', D.srate, 'bla');
nC = size(DS,2);
LN = 1000000;
Tfbotm = size(T,2)/nC;
Cest = mysort.noise.Covest2(DS, 'maxLag', Tfbotm, 'maxSamples', LN, 'noiseEpochs', noiseEpochs);

% Cest.CCol(1:4,:) = Cest.CCol(1:4,:) + diag(diag(Cest.CCol(1:4,:)));
% Cest.CCol = Cest.CCol/2;
matchedUnit = 4;
cutLeft = 18;
% Cest =();
tT =  mysort.wf.v2t(T, nC);
botm = mysort.sorters.BOTM(Cest, tT, 'adaptOnInit', 1);
botm.samplesPerSecond = 20000;
botm.DH = DS;    
GTSSC = mysort.spiketrain.SpikeSortingContainer('GTIntra', [matchedUnit*ones(length(D.gt_spiketrain), 1) D.gt_spiketrain(:)], 'templateWfs', tT(:,:,matchedUnit), 'templateCutLeft', cutLeft);
DS.addSpikeSorting(GTSSC);
SDA = mysort.plot.SliderDataAxes({DSI DS botm}, 'channelSpacers', [0 1200  0], 'timeIn', 'sec');
SDA.setXLim((3.22223e+06 + [0 3200])/20000);
SDA.update();

% SDA.bAvoidReentrant = 1;
% SDA.bIsVisible = 0;
% mysort.plot.savefig(gcf, 'HarrisBOTM_Export_forPaperFigure', 'fig', 0, 'ai',1)
return
% Good spots:
% 71500, L 8000
% 133500, L 2000
% 22500, L 64000
% 3.22223e+06, L 3200

%%
figure
ah = subplot(3,1,1)
hold on
ah(2) = subplot(3,1,2)
hold on
ah(3) = subplot(3,1,3)
hold on
for  i=1:nT
%     plot(edges, nEucl, 'color', mysort.plot.vectorColor(i), 'linewidth', 2)
    plot(ah(1),edges, nEuclCor(:,i), '-', 'color', mysort.plot.vectorColor(i), 'linewidth', 3);
    plot(ah(2),edges, nConvCor(:,i), '-', 'color', mysort.plot.vectorColor(i), 'linewidth', 3);
    plot(ah(3),botm_edges, nBotmCor(:,i), '-', 'color', mysort.plot.vectorColor(i), 'linewidth', 3);
end
linkaxes(ah, 'xy');
set(ah(3), 'xlim', [0 10e6]);
plot(0,0,'.', 'color', mysort.plot.vectorColor(correctId), 'markersize', 20)

%%
scalefactor = 1/mean(sqrt(diag(C))); 
scalefactor = 1/max((sqrt(diag(C))));
mysort.plot.figure('w', 1000, 'h', 1000);
p=0; nR = 5; nC=2;
% Eucl
p=p+1; ah(p) = subplot(nR,nC,p);
dims = size(T,2);
width = edges(2)-edges(1);
scalefactor = 1/mean(sqrt(diag(C_)));
% Eucl
% plot(ah(p), edges, nEucl(:,correctIdx)/(sum(nEucl(:,correctIdx))*width), 'b', 'linewidth', 2);
hold on
plot(ah(p), edges, nEuclCor(:,correctIdx)/(sum(nEuclCor(:,correctIdx))*width), 'g', 'linewidth', 2);
e = norm(T(correctIdx,:));
% theo = 10e4*scalefactor^2*ncx2pdf(scalefactor^2*(edges), Tf, (scalefactor*e)^2);
theo_correct = scalefactor^2*chi2pdf(scalefactor^2*(edges/10), dims);
theo_noise   = scalefactor^2*ncx2pdf(scalefactor^2*edges/10, dims, (scalefactor*e)^2);
plot(ah(p), edges/10, theo_correct, 'g:', 'linewidth', 2);
plot(ah(p), edges/10, theo_noise, 'r:', 'linewidth', 2);
plot(ah(p), edges, nEuclNoise(:,correctIdx)/(sum(nEuclNoise(:,correctIdx))*width), 'r', 'linewidth', 2);
E = @(thr) sum(nEuclCor(:,correctIdx)<thr) 
p=p+1; ah(p) = subplot(nR,nC,p);
[mx mxidx] = min(EuclCor,[],2);
n = histc(mxidx, [0:1:12]);
classPerfEucl = sum(mxidx==correctIdx)/length(mxidx);
bar(1:12, n(1:12))
title('Eucl');

% Maha
dims = 139;
p=p+1; ah(p) = subplot(nR,nC,p);
width = botm_edges(2)-botm_edges(1);
% plot(ah(p), botm_edges, nMaha(:,correctIdx)/(sum(nMaha(:,correctIdx))*width), 'b', 'linewidth', 2);
hold on
plot(ah(p), botm_edges, nMahaCor(:,correctIdx)/(sum(nMahaCor(:,correctIdx))*width), 'g', 'linewidth', 2);
plot(ah(p), botm_edges, nMahaNoise(:,correctIdx)/(sum(nMahaNoise(:,correctIdx))*width), 'r', 'linewidth', 2);
e = sqrt((T(correctIdx,:)/C_)*T(correctIdx,:)');
scalefactor = 1;
% theo = 10e4*scalefactor^2*ncx2pdf(scalefactor^2*(edges), Tf, (scalefactor*e)^2);
theo_correct = scalefactor^2*chi2pdf(scalefactor^2*botm_edges, dims);
theo_noise   = scalefactor^2*ncx2pdf(scalefactor^2*botm_edges, dims, (scalefactor*e)^2);
plot(ah(p), botm_edges, theo_correct, 'g:', 'linewidth', 2);
plot(ah(p), botm_edges, theo_noise, 'r:', 'linewidth', 2);

p=p+1; ah(p) = subplot(nR,nC,p);
[mx mxidx] = min(MahaCor,[],2);
n = histc(mxidx, [0:1:12]);
classPerfMaha = sum(mxidx==correctIdx)/length(mxidx);
bar(1:12, n(1:12))
title('Maha');

% Conv
p=p+1; ah(p) = subplot(nR,nC,p);
% plot(ah(p), edges, nConv(:,correctIdx)/(sum(nConv(:,correctIdx))*width), 'b', 'linewidth', 2);
hold on
plot(ah(p), edges, nConvCor(:,correctIdx)/(sum(nConvCor(:,correctIdx))*width), 'g', 'linewidth', 2);
plot(ah(p), edges, nConvNoise(:,correctIdx)/(sum(nConvNoise(:,correctIdx))*width), 'r', 'linewidth', 2);

p=p+1; ah(p) = subplot(nR,nC,p);
[mx mxidx] = max(ConvCor,[],2);
classPerfConv = sum(mxidx==correctIdx)/length(mxidx);
n = histc(mxidx, [0:1:12]);
bar(1:12, n(1:12))
title('Conv');

% Matched
p=p+1; ah(p) = subplot(nR,nC,p);
% plot(ah(p), botm_edges, nMatched(:,correctIdx)/sum(nMatched(:,correctIdx)), 'b', 'linewidth', 2);
hold on
plot(ah(p), botm_edges, nMatchedCor(:,correctIdx)/sum(nMatchedCor(:,correctIdx)), 'g', 'linewidth', 2);
plot(ah(p), botm_edges, nMatchedNoise(:,correctIdx)/sum(nMatchedNoise(:,correctIdx)), 'r', 'linewidth', 2);

p=p+1; ah(p) = subplot(nR,nC,p);
[mx mxidx] = max(MatchedCor,[],2);
classPerfMatched = sum(mxidx==correctIdx)/length(mxidx);
n = histc(mxidx, [0:1:12]);
bar(1:12, n(1:12))
title('Matched');

% BOTM
p=p+1; ah(p) = subplot(nR,nC,p);
% plot(ah(p), botm_edges, nBotm(:,correctIdx)/sum(nBotm(:,(:,correctIdx)), 'b', 'linewidth', 2);
hold on
plot(ah(p), botm_edges, nBotmCor(:,correctIdx)/sum(nBotmCor(:,correctIdx)), 'g', 'linewidth', 2);
plot(ah(p), botm_edges, nBotmNoise(:,correctIdx)/sum(nBotmNoise(:,correctIdx)), 'r', 'linewidth', 2);

p=p+1; ah(p) = subplot(nR,nC,p);
[mx mxidx] = max(botmCor,[],2);
classPerfBotm = sum(mxidx==correctIdx)/length(mxidx);
n = histc(mxidx, [0:1:12]);
bar(1:12, n(1:12))
title('Botm');

xlabel('Assigned Class', 'fontsize', 14)
set(ah, 'fontsize', 14);
% linkaxes(ah([1 3 5]), 'xy');
set(ah([1 5]), 'xlim', [0 10e6]);
set(ah([3]), 'xlim', [0 500]);
set(ah([7 9]), 'xlim', [-250 250]);

%%
figure
plot(S.spikeCut.wfs', 'k')
hold on
plot(T(1,:), 'g')

%%
figure
plot(correctSpikes', 'k')
hold on
plot(T(5,:), 'g')

%%
figure
bar([classPerfEucl classPerfMaha classPerfConv classPerfMatched classPerfBotm])


%%
dims = size(T,2);
width = edges(2)-edges(1);
scalefactor = 1/mean(sqrt(diag(C_)));
mysort.plot.figure('w', 1000, 'h', 1000);
% Eucl
ah = axes;
p = 1;
% plot(ah(p), edges, nEucl(:,correctIdx)/(sum(nEucl(:,correctIdx))*width), 'b', 'linewidth', 2);
hold on
plot(ah(p), edges, nEuclCor(:,correctIdx)/(sum(nEuclCor(:,correctIdx))*width), 'g', 'linewidth', 2);
e = norm(T(correctIdx,:));
% theo = 10e4*scalefactor^2*ncx2pdf(scalefactor^2*(edges), Tf, (scalefactor*e)^2);
theo_correct = scalefactor^2*chi2pdf(scalefactor^2*(edges), dims);
theo_noise   = scalefactor^2*ncx2pdf(scalefactor^2*edges, dims, (scalefactor*e)^2);
plot(ah(p), edges, theo_correct, 'g:', 'linewidth', 2);
plot(ah(p), edges, theo_noise, 'r:', 'linewidth', 2);
plot(ah(p), edges, nEuclNoise(:,correctIdx)/(sum(nEuclNoise(:,correctIdx))*width), 'r', 'linewidth', 2);
set(gca, 'xlim', [0 1.2*10e6])