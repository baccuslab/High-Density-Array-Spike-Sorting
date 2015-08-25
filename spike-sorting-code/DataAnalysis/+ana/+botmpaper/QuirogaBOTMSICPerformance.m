benchmarks = {'Easy1', 'Easy2', 'Difficult1', 'Difficult2'};
noiselevels = {'005', '01', '015', '02'};
initTf = 81;
initCutLeft = -5;
Tf = 71;
cutLeft = 15;
dpath = ana.botmpaper.E10_quirogaDataPath();
samplesPerSecond = 32000;
Eucl  = @ana.botmpaper.Eucl;
Maha  = @ana.botmpaper.Maha;
Conv  = @ana.botmpaper.Conv;
Match = @ana.botmpaper.Match;
Botm  = @ana.botmpaper.Botm;
global savefig__;
savefig__ = 1;
figSavePath = 'C:\Dropbox\PaperBOTM1\DetectorFigures\';
cols = {'b', 'g', 'r'};
pd = pdefs();
runName = 'V03_DL.5';
runName = 'V02';
runName = 'V04_BOTMwoSIC';
runName = 'V04_BOTMwoSIC_DL';
outPath = fullfile(pd.serverData, 'Quiroga', 'BOTMPaper1SICRun');


%%
bUseLambdaDL = 0; % This switched between DL with condition and DL with lambda
lambda = .3;
resultT = zeros(4*4, 3);
resultT0 = zeros(4*4, 3);
resultT_wc = zeros(4*4, 3);
resultT0_wc = zeros(4*4, 3);
SPIKEPRIOR = .001;
NOISEPRIOR = 1-3*SPIKEPRIOR;
refactory = 10;
refactoryPeaks = 5;
upsample = 4;
maxIter = 6;

%% DO THE OLD ANALYSIS WITH THE GT DERIVED TEMPLATES
if 0
    outPathRun = fullfile(outPath, runName);
    if ~exist(outPath, 'file'); mkdir(outPath); end
    if ~exist(outPathRun, 'file'); mkdir(outPathRun); end
    bUSE_WAVE_CLUS_TEMPALTES = false;
    filei = 0;
    for b = 1:length(benchmarks)
        for n = 1:length(noiselevels)
            filei = filei+1;
            fname = fullfile(outPathRun, sprintf('OvpRunR_%d_%d.mat', b, n));
            fnameY = fullfile(outPathRun, sprintf('OvpRunR_%d_%d_Y.mat', b, n));
            if 1%~exist(fname, 'file')
                ana.botmpaper.QuirogaBOTMSICPerformanceBody()
                save(fname, 'R_sic', 'gdf_sic', 'gdf_sic_up', 'lambda', 'upsample', 'Mup', 'M', 'refactory', 'maxIter', 'bUseLambdaDL', 'tgdfshifted')
    %             save(fnameY, 'Ysic');
                resultT0(filei,:) = [sum(R_sic{2}.detErr) sum(R_sic{2}.nCL) sum(R_sic{2}.totErr)];
                resultT(filei,:) = [sum(R_sic{end}.detErr) sum(R_sic{end}.nCL) sum(R_sic{end}.totErr)];
            end
        end
    end
    resultT
    sum(resultT)
end

%% DO THE NEW ANALYSIS WITH THE WAVE_CLUS DERIVED TEMPLATES
bUSE_WAVE_CLUS_TEMPALTES = true;
outPathRun = fullfile(outPath, [runName '_Wave_CLUS_TEMPALTES']);
if ~exist(outPath, 'file'); mkdir(outPath); end
if ~exist(outPathRun, 'file'); mkdir(outPathRun); end
filei = 0;

for b = 1:length(benchmarks)
    for n = 1:length(noiselevels)
        filei = filei+1;
        fname = fullfile(outPathRun, sprintf('OvpRunR_%d_%d.mat', b, n));
        fnameY = fullfile(outPathRun, sprintf('OvpRunR_%d_%d_Y.mat', b, n));
        if 1%~exist(fname, 'file')
            ana.botmpaper.QuirogaBOTMSICPerformanceBody()
%             save(fname, 'R_sic', 'gdf_sic', 'gdf_sic_up', 'lambda', 'upsample', 'Mup', 'M', 'refactory', 'maxIter', 'bUseLambdaDL', 'tgdfshifted')
%             save(fnameY, 'Ysic');
            resultT0_wc(filei,:) = [sum(R_sic{2}.detErr) sum(R_sic{2}.nCL) sum(R_sic{2}.totErr)];
            resultT_wc(filei,:) = [sum(R_sic{end}.detErr) sum(R_sic{end}.nCL) sum(R_sic{end}.totErr)];
        end
    end
end
resultT_wc
sum(resultT_wc)

%%
NSpikes = [];
i = 1;
for b=1:length(benchmarks)
    for n=1:length(noiselevels)
        fprintf('\n###################\nProcessing %d - %d \n', b, n)
        quirogaFile = sprintf('C_%s_noise%s.mat', benchmarks{b}, noiselevels{n});
        data = load(fullfile(simulator_path, quirogaFile));
        NSpikes(i,1) = sum(cellfun(@(x) length(x), data.spike_class));
        i = i+1;
    end
end
TP = NSpikes-resultT_wc(:,1);
FN = resultT_wc(:,1);
DETPerf =  100*TP./(TP+FN);
TP = NSpikes - resultT_wc(:,2);
CN = resultT_wc(:,2);
CLAPerf = 100*TP./(TP+CN);
TotalPerf = (DETPerf + CLAPerf)/2;
round(10*([DETPerf CLAPerf TotalPerf]))
% resultT_wc =
% 
%           11           2          13
%            4           2           6
%            8           0           8
%            9           3          12
%            2           2           4
%            6           2           8
%            4           6          10
%            6           7          13
%            2          18          20
%           18          12          30
%            9          17          26
%           20          10          30
%            8           9          17
%            5           9          14
%            8           9          17
%           36        1124        1160
return

%% COMPILE RESULTS
resultT = zeros(4*4, 3);
resultT0 = zeros(4*4, 3);
i = 0;
for b = 1:length(benchmarks)
    for n = 1:length(noiselevels)
        i=i+1;
        fname = fullfile(outPathRun, sprintf('OvpRunR_%d_%d.mat', b, n));
%         try
            clear R_sic;
            load(fname);
            resultT(i,:) = [sum(R_sic{end}.detErr) sum(R_sic{end}.nCL) sum(R_sic{end}.totErr)];
            resultT0(i,:) = [sum(R_sic{2}.detErr) sum(R_sic{2}.nCL) sum(R_sic{2}.totErr)];
%         end
    end
end    
resultT

return



%%
b = 2
n = 4
fname = fullfile(outPathRun, sprintf('OvpRunR_%d_%d.mat', b, n));
load(fname);
mysort.plot.printEvaluationTable(R_sic{5})

%% PLOT FNs OVERLAPS !
def = mysort.util.defs();
stidx = 1;
fpidx = find(R_sic{4}.spikeLabel1{stidx} == def.fno);
fp_times = R_sic{4}.St1{stidx}(fpidx);
W = 1000;
for k=1:min(6,length(fp_times))
    fpt = fp_times(k);
    xrange = (fpt-W):fpt+W;
    xrange_up = upsample*(fpt-W):upsample*(fpt+W);
    tgdfidx =  (fpt-W) < tgdfshifted(:,2)   &    tgdfshifted(:,2) < (fpt+W);
   
    gdf3idx =  (fpt-W) < gdf_sic{3}(:,2)   &    gdf_sic{3}(:,2) < (fpt+W);
    
    mysort.plot.figure([800, 1000]);
    nR = 4; nC=1; p=0; ah = [];
    p=p+1; ah(p) = subplot(nR, nC, p);
    plot(xrange, GT.X(xrange)')
    hold on
    mysort.plot.spiketrain(tgdfshifted(tgdfidx,:), 'yoffset', 0, 'colorList', cols)
    axis tight
    
    p=p+1; ah(p) = subplot(nR, nC, p);
    plot(xrange_up/upsample, Yup(:,xrange_up)')
%     hold on
%     mysort.plot.spiketrain(tgdfshifted(tgdfidx,:), 'yoffset', 0, 'colorList', cols)
%     axis tight

    for i=2:3
        p=p+1; ah(p) = subplot(nR, nC, p);
        gdf_idx =  (fpt-W) < gdf_sic{i}(:,2)   &    gdf_sic{i}(:,2) < (fpt+W);
        plot(xrange_up/upsample, Ysic{i}(:,xrange_up)')
        hold on
        mysort.plot.spiketrain(gdf_sic{i}(gdf_idx,:), 'yoffset', 0, 'colorList', cols)
        axis tight
    end

    linkaxes(ah, 'x');
    mysort.plot.dataCursorMode();
end

%% PLOT FPs!
def = mysort.util.defs();
stidx = 3;
fpidx = find(R_sic{4}.spikeLabel2{stidx} == def.fp);
fp_times = R_sic{4}.St2{stidx}(fpidx);
W = 1000;
for k=1:min(1,length(fp_times))
    fpt = fp_times(k);
    xrange = max(1,(fpt-W)):fpt+W;
    xrange_up = upsample*max(1,(fpt-W)):upsample*(fpt+W);
    tgdfidx =  (fpt-W) < tgdfshifted(:,2)   &    tgdfshifted(:,2) < (fpt+W);
   
    gdf3idx =  (fpt-W) < gdf_sic{3}(:,2)   &    gdf_sic{3}(:,2) < (fpt+W);
    
    mysort.plot.figure([800, 1000]);
    nR = 6; nC=1; p=0; ah = [];
    p=p+1; ah(p) = subplot(nR, nC, p);
    plot(xrange, GT.X(xrange)')
    hold on
    mysort.plot.spiketrain(tgdfshifted(tgdfidx,:), 'yoffset', 0, 'colorList', cols)
    axis tight
    
    p=p+1; ah(p) = subplot(nR, nC, p);
    plot(xrange_up/upsample, Ysic{1}(:,xrange_up)')
%     hold on
%     mysort.plot.spiketrain(tgdfshifted(tgdfidx,:), 'yoffset', 0, 'colorList', cols)
%     axis tight

    for i=2:5
        p=p+1; ah(p) = subplot(nR, nC, p);
        gdf_idx =  (fpt-W) < gdf_sic{i}(:,2)   &    gdf_sic{i}(:,2) < (fpt+W);
        plot(xrange_up/upsample, Ysic{i}(:,xrange_up)')
        hold on
        mysort.plot.spiketrain(gdf_sic{i}(gdf_idx,:), 'yoffset', 0, 'colorList', cols)
        axis tight
    end

    linkaxes(ah, 'x');
    mysort.plot.dataCursorMode();
end
%%
mysort.plot.alignment(R)
%% PLOT CLASSIFICATION ERRORS, Overlapping SPIKES!
def = mysort.util.defs();
stidx = 3;
fpidx = find(R_sic{4}.spikeLabel2{stidx} == def.clo);
fp_times = R_sic{4}.St2{stidx}(fpidx);
W = 1000;
for k=1:min(6,length(fp_times))
    fpt = fp_times(k);
    xrange = max(1,(fpt-W)):fpt+W;
    xrange_up = upsample*max(1,(fpt-W)):upsample*(fpt+W);
    tgdfidx =  (fpt-W) < tgdfshifted(:,2)   &    tgdfshifted(:,2) < (fpt+W);
   
    gdf3idx =  (fpt-W) < gdf_sic{3}(:,2)   &    gdf_sic{3}(:,2) < (fpt+W);
    
    mysort.plot.figure([800, 1000]);
    nR = 4; nC=1; p=0; ah = [];
    p=p+1; ah(p) = subplot(nR, nC, p);
    plot(xrange, GT.X(xrange)')
    hold on
    mysort.plot.spiketrain(tgdfshifted(tgdfidx,:), 'yoffset', 0, 'colorList', cols)
    axis tight
    
    p=p+1; ah(p) = subplot(nR, nC, p);
    plot(xrange_up/upsample, Yup(:,xrange_up)')
%     hold on
%     mysort.plot.spiketrain(tgdfshifted(tgdfidx,:), 'yoffset', 0, 'colorList', cols)
%     axis tight

    for i=2:3
        p=p+1; ah(p) = subplot(nR, nC, p);
        gdf_idx =  (fpt-W) < gdf_sic{i}(:,2)   &    gdf_sic{i}(:,2) < (fpt+W);
        plot(xrange_up/upsample, Ysic{i}(:,xrange_up)')
        hold on
        mysort.plot.spiketrain(gdf_sic{i}(gdf_idx,:), 'yoffset', 0, 'colorList', cols)
        axis tight
    end

    linkaxes(ah, 'x');
    mysort.plot.dataCursorMode();
end

