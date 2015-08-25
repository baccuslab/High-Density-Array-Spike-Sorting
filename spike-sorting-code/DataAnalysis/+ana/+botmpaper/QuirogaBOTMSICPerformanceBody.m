savePath = fullfile(dpath, benchmarks{b}, noiselevels{n});
if ~exist(savePath, 'file')
    mkdir(savePath);
end      

pic = 0;
fprintf('Starting with %s %s\n', benchmarks{b}, noiselevels{n});
close all
quirogaFile = sprintf('C_%s_noise%s.mat', benchmarks{b}, noiselevels{n});        
GT = ana.botmpaper.E10_preprocessing(dpath, quirogaFile, initCutLeft, initTf, dpath);

bufferFile = fullfile(savePath, 'botmOutputDistributionsThrehold.mat');
        
% RUN QuirogaOptimalThreshold first to save all buffered files.
load(bufferFile);
tgdfshifted = tgdf;
tgdfshifted(:,2) = tgdfshifted(:,2)+cutLeft;
% T(4,:) = mysort.util.shiftSubtract(T(1,:), -T(2,:), 5);

%% COMPUTE FILTER OUTPUTS
if bUSE_WAVE_CLUS_TEMPALTES
    Tf_cut = 81;
    simulator_path = ana.botmpaper.E10_quirogaDataPath();
    wc_file = sprintf('times_C_%s_noise%s.mat', benchmarks{b}, noiselevels{n});
    sorted = load(fullfile(simulator_path, wc_file));
    sorted_gdf = sorted.cluster_class;
    sorted_gdf(sorted_gdf(:,1) == 0,:) = [];
    
    epochs = [sorted_gdf(:,2) sorted_gdf(:,2)+Tf_cut-1]-cutLeft;
    spikes = mysort.epoch.extractWaveform(GT.X, epochs);
    T = mysort.util.calculateClassMeans(spikes, sorted_gdf(:,1));
else
    T = GT.templates;
end
C = GT.C;
nT = size(T,1);
nC = 1;
Y = zeros(size(GT.X));
Tf2 = floor(size(T,2)/2);
M = zeros(2*Tf2+1, nT, nT);
if bUseLambdaDL
    C_ = (1-lambda)*C + lambda*diag(diag(C));
else
    C_ = ana.botmpaper.diagonalLoading(GT.C, 10000);
end
for i=1:nT
    t = T(i,:);
                    t_ = mysort.util.embedTime2embedChan(t, nC);
                    CCol = mysort.noise.Cte2Ccol(C_, nC);
                    f_ = matlabfilecentral.block_levinson(t_(:), CCol)';
                    f = mysort.util.embedChan2embedTime(f_, nC);
    Y(i,:)  = conv(GT.X, flipud(f(:)), 'same') - .5*t*f' + log(SPIKEPRIOR);
    for j=1:nT
        M(:,i,j) = conv(T(j,:), flipud(f(:)), 'same');
    end
end

%% PLOT XI vs F
if 0
figure;
Tf2 = floor(size(M,1)/2);
range = -Tf2:Tf2;
for f1 = 1:nT
    for f2 = 1:nT
        pIdx = (f1-1)*nT + f2;
        subplot(nT,nT,pIdx);
        x = squeeze(M(:,f1,:));
        plot(range,x, '-.', 'linewidth', 2);
        hold on
        plot([0 0], [min(x(:)) max(x(:))], 'k:');
    end
end
end
%%
Yup = resample(Y', upsample, 1)';
Mup = mysort.wf.tResample(M, upsample, 1);
refactoryup = refactory*upsample;
refactoryPeaksUp = refactoryPeaks*upsample;
Ysic = {}; Ysic{1} = Yup;
gdf_sic_up = {};
gdf_sic = {}; gdf_sic{1} = [];
R_sic = {};
for it = 2:maxIter
    [gdf_sic_up{it} Ysic{it}] = ana.botmpaper.QuirogaBOTMSICPerformanceDOSIC(Ysic{it-1}, Mup, refactoryup, refactoryPeaksUp, NOISEPRIOR, SPIKEPRIOR);
    gdf = gdf_sic_up{it}; gdf(:,2) = round(gdf(:,2)/upsample);
    gdf_sic{it} = sortrows([gdf_sic{it-1}; gdf],2);
    R_sic{it} = mysort.spiketrain.alignGDFs(tgdfshifted, gdf_sic{it}, 25, 20, 25);
    mysort.plot.printEvaluationTable(R_sic{it}, 'ignoreEmptySorted', false)
end




