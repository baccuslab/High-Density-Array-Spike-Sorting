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
nT = size(T,1);
Y = zeros(size(GT.X));
Tf2 = floor(size(T,2)/2);
M = zeros(2*Tf2+1, nT, nT);
C_ = (1-lambda)*C + lambda*diag(diag(C));
for i=1:nT
    t = T(i,:);
    f = t/C_;
    Y(i,:)  = conv(GT.X, flipud(f(:)), 'same') - .5*t*f' + log(.001);
    for j=1:3
        M(:,i,j) = conv(T(j,:), flipud(f(:)), 'same');
    end
end
% maxTau = 5;
% DC = mysort.sorters.DiscriminantFunctionContainer(Y', M, maxTau);
%%
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

%%

[DMax Ids OvpIndex] = mysort.sorters.DiscriminantFunctionMaximums(Y', M, 25, resampleP);

%%
refactory = 25;
[gdf_resolved gdf_raw] = mysort.sorters.DiscriminantFunctions2Gdf(Y', DMax, Ids, OvpIndex, refactory, 10);
% gdf_raw(:,2) = gdf_raw(:,2)-cutLeft;
% gdf_resolved(:,2) = gdf_resolved(:,2)-cutLeft;

% figure
% plot(T([1:2],:)')
% hold on
% plot()
% 
% figure; hist(gdf_raw(:,1), 120)

gdf_raw_ = gdf_raw;
gdf_raw_(gdf_raw(:,1)>4,1) = 4;
R_ovp_ = mysort.spiketrain.alignGDFs(tgdfshifted, gdf_raw_, 10, 3, 25);
mysort.plot.printEvaluationTable(R_ovp_)

gdf_ovp = gdf_resolved;
R_ovp = mysort.spiketrain.alignGDFs(tgdfshifted, gdf_ovp, 10, 3, 25);
mysort.plot.printEvaluationTable(R_ovp)