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
% Good spots:
% 71500, L 8000
% 133500, L 2000
% 22500, L 64000
% 3.22223e+06, L 3200
