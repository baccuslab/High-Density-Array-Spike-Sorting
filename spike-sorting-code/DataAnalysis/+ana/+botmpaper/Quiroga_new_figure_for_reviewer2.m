benchmarks = {'Easy1', 'Easy2', 'Difficult1', 'Difficult2'};
noiselevels = {'005', '01', '015', '02'};
initTf = 81;
initCutLeft = -5;
Tf = 61;
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
%%
b = 4
n = 4
close all
savePath = fullfile(dpath, benchmarks{b}, noiselevels{n});
if ~exist(savePath, 'file')
    mkdir(savePath);
end        
pic = 0;
fprintf('Starting with %s %s\n', benchmarks{b}, noiselevels{n});
close all
quirogaFile = sprintf('C_%s_noise%s.mat', benchmarks{b}, noiselevels{n});
GT = ana.botmpaper.E10_preprocessing(dpath, quirogaFile, initCutLeft, initTf, dpath);
gt_spikes = mysort.epoch.extractWaveform(GT.X, [GT.gdf(:,2)-initCutLeft GT.gdf(:,2)-initCutLeft+initTf-1]);
T = mysort.util.calculateClassMeans(gt_spikes, GT.gdf(:,1));

[alignedTemplates tau] = mysort.wf.vAlignOnMax(T, 1, 'truncate', 0);
alignedTemplates = alignedTemplates(:, 5:75);
[maxVal maxIdx] = max(alignedTemplates(1,:));
tgdf = GT.gdf;
for i=1:3
    tgdf(GT.gdf(:,1)==i,2) = GT.gdf(GT.gdf(:,1)==i,2) + maxIdx - tau(i);
end
gt_spikes = mysort.epoch.extractWaveform(GT.X, [tgdf(:,2)-cutLeft tgdf(:,2)-cutLeft+Tf-1]);
T = mysort.util.calculateClassMeans(gt_spikes, tgdf(:,1));

% Figure 1
mysort.plot.figure('w',1000,'h',1000); 
subplot(3,1,1)
plot(GT.templates');
title('original Templates, cut');
subplot(3,1,2)
plot(alignedTemplates');
title('aligned Templates, cut');
subplot(3,1,3)
plot(gt_spikes');
hold on
plot(T', 'linewidth', 3, 'color', 'k');
title('realigned gt Spikes');        
pic = pic+1; mysort.plot.savefig(gcf, fullfile(savePath, [sprintf('%03d', pic) 'Templates']), 'fig', 0);

% Noise estimation
spikeEpochs = mysort.epoch.merge([tgdf(:,2)-20 tgdf(:,2)+80]);
noiseEpochs = mysort.epoch.flip(spikeEpochs, size(GT.X,2));
noiseEpochs = mysort.epoch.removeShort(noiseEpochs, Tf);
DS = mysort.ds.Matrix(GT.X', 32000, 'Q', 1, 1);
noiseSnippets = DS.getWaveform(noiseEpochs(:,1), 0, Tf);
nNoise = size(noiseSnippets,1);

XCC = mysort.noise.XCorrContainer(DS, Tf-1, 'noiseEpochs', noiseEpochs);
Cest = XCC.getCte4Channels(1:size(DS,2));
dCest = diag(diag(Cest));
C = cov(noiseSnippets);
dC = diag(diag(C));
targetC = dC;
targetCest = dCest;        
% see Pope 2008, "Shrinkage estimation of the power spectrum covariance matrix", eq.10
meanNoise = mean(noiseSnippets,1);
meanFreeNoise = noiseSnippets - repmat(meanNoise,nNoise,1);
varC  = zeros(size(C));
covTC = varC;
for k=1:nNoise
    Wijk = meanFreeNoise(k,:)'*meanFreeNoise(k,:);
    Tijk = diag(diag(Wijk));
    varC = varC + (Wijk - C).^2;
    covTC = covTC + (Wijk - C).*(Tijk - targetC);
end
varC = varC * nNoise/(nNoise-1)^3;
covTC = covTC * nNoise/(nNoise-1)^3;    
lambdaopt = sum(varC(:) - covTC(:))   /  sum(sum((C-targetC)).^2);
lambda = .5;
C_ = (1-lambda)*C + lambda*targetC;
C_opt = (1-lambdaopt)*C + lambdaopt*targetC;
Cest_ = (1-lambda)*Cest + lambda*targetCest; 
CDL = ana.botmpaper.diagonalLoading(C, 10000);
CestDL = ana.botmpaper.diagonalLoading(Cest, 10000);

%
mysort.plot.figure('w',1000,'h',1000);
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
pic=pic+1; mysort.plot.savefig(gcf, fullfile(savePath, [sprintf('%03d_', pic) 'covs']));   


% Projection estimation
H = cov(gt_spikes);
[VH DH] = eig(H);
[VC DC] = eig(C);
[VC_ DC_] = eig(C_);
dH = diag(DH);
largeEV = dH > .001*max(dH);
Proj = VH(:, largeEV);
fprintf('Chosen number of dimensions: %d\n', size(Proj,2));   


CP = cov(noiseSnippets*Proj);
CPest = Proj'*Cest*Proj;
HP = cov(gt_spikes*Proj);
TP = T*Proj;      

%
figure
plot(diag(DH))
hold on
plot(diag(DC), 'r')
plot(diag(DC_), 'r:')
legend('H', 'C', 'C.')
title('Eigenvalues of covariance matrices');
pic=pic+1; mysort.plot.savefig(gcf, fullfile(savePath, [sprintf('%03d_', pic) 'eigenvalues_of_covs']));        

%
DATASETS = {gt_spikes, T, C, 'corr, C', 0, []
            noiseSnippets,  T, C, 'N, C', 0, []
            gt_spikes, T, dC, 'corr, dC', 0, []
            noiseSnippets,  T, dC, 'N, dC', 0, []
            gt_spikes, T, C_opt, 'corr, C_opt', 0, []
            noiseSnippets,  T, C_opt, 'N, C_opt', 0, []
            gt_spikes, T, C_, 'corr, C_', 0, []
            noiseSnippets,  T, C_, 'N, C_', 0 , []               
            gt_spikes, TP, CP, 'corr, CP', 0, Proj
            noiseSnippets, TP, CP, 'N, CP', 0 , Proj    

            gt_spikes, T, Cest, 'corr, Cest', 0 , [] 
            noiseSnippets,  T, Cest, 'N, Cest', 0 , [] 
            gt_spikes, T, dCest, 'corr, dC', 0 , [] 
            noiseSnippets,  T, dCest, 'N, dC', 0 , []   
            gt_spikes, T, Cest_, 'corr, Cest_', 0 , [] 
            noiseSnippets,  T, Cest_, 'N, Cest_', 0 , [] 
            gt_spikes, TP, CPest, 'corr, CPest', 0, Proj
            noiseSnippets, TP, CPest, 'N, CPest', 0, Proj
            gt_spikes, T, CDL, 'corr, CDL', 0 , [] 
            noiseSnippets,  T, CDL, 'N, CDL', 0 , []
            gt_spikes,  T, CestDL, 'corr, CestDL', 0, []
            noiseSnippets,  T, CestDL, 'N, CestDL', 0, []};

datasetNames = DATASETS(:,4);

methods = { Eucl, 'noC', 'minDetection', 'Eucl'
            Maha, 'C',   'minDetection', 'Maha'
            Conv, 'noC', 'maxDetection', 'Conv'
            Match, 'C',  'maxDetection', 'Match'
            Botm, 'C',   'maxDetection', 'Botm'
            Botm, 'C',   'maxDetection', 'BotmOpt'};

% Compute Classifications
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
            else
                RES(d).TM(:,m,i) = fun(DATASETS{d,1}, lT(i,:), DATASETS{d,6});
                RES(d).cputime(m,i) = toc(t1)+DATASETS{d,5};
            end
        end
        fprintf('Estimating performance...\n');
        if strcmp(methods{m,3}, 'minDetection')
            [mx mxidx] = min(squeeze(RES(d).TM(:,m,:)),[],2);
        else
            [mx mxidx] = max(squeeze(RES(d).TM(:,m,:)),[],2);
        end            
        RES(d).classifications(:,m) = histc(mxidx, [0.5:1:(nT+0.5)]);
        if ~isempty(strfind(DATASETS{d,4}, 'corr'))
            RES(d).performance(m) = sum(mxidx==tgdf(:,1))/length(mxidx);   
        else
            RES(d).performance(m) = 0;   
        end

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


        
  

%%
dataSetIdxCorr = 11;
dataSetIdxN = 12;

dataSetIdxCorr = 15;
dataSetIdxN = 16;

DATASETS{dataSetIdxCorr, 4}
DATASETS{dataSetIdxN, 4}
methods
R_CORR = RES(dataSetIdxCorr);
R_N    = RES(dataSetIdxN);
IDs = GT.gdf(:,1);

botmI = 5;
D_matchedI = 4;

TMbotm = squeeze(R_CORR.TM(:,botmI,:));
TMmatc = squeeze(R_CORR.TM(:,D_matchedI,:));
MI = min(TMbotm(:));
MA = max(TMbotm(:));



%%
% X = TMbotm;
% MI = min(X);
% MA = max(X);
% 
% 
% 
% % for i=1:3
% %     idx1 = IDs==i; %ismember(IDs, [1 2]);
% % %     M(i,:) = mean(TMbotm(idx1,:));
% % end
% % D1 = M(1,:) - M(2,:);
% % D1 = D1/norm(D1);
% % D2 = M(2,:) - M(3,:);
% % 
% % D2 = D2 - D1*(D2*D1');
% % D2 = D2/norm(D2);
% % 
% % XP = X*[D1' D2'];
% 
% mysort.plot.figure([400 400]);
% nR = 1; nC = 1; 
% cols = {'b', 'g', 'r'};
% ahx = mysort.plot.subplots([nR nC], 'matrix', 1, 'spacerX', .1);
% for i=1:3
%     idx1 = IDs==i; %ismember(IDs, [1 2]);
%     plot3(ahx,X(idx1,1), X(idx1,2), X(idx1,3),'.', 'color', cols{i});
% end
% 
% return


%%
XX = {TMbotm TMmatc};

XLIMS = {[-20 70] [-20 70] [-25 75]};
YLIMS = {[-20 70] [-20 70] [-20 70]};
     
     
BotmConsts = TMbotm(1,:) - TMmatc(1,:);
Pairs = [1 2
         1 3
         2 3];
for k=1:1 %length(XX)
    X = XX{k};
    MI = min(X);
    MA = max(X);

    mysort.plot.figure([1000 270]);
    nR = 1; nC = 3; 
    ahx = mysort.plot.subplots([nR nC], 'matrix', 1, 'spacerX', .1);
    p = 1;
    for ii=1:3
        i = Pairs(ii,1);
        kk = Pairs(ii,2);
        idx1 = IDs==i; %ismember(IDs, [1 2]);
        idx2 = IDs == kk;
        idx3 = ~idx1 & ~idx2;
            
        mii= min(MI([i kk]));
        maa= max(MA([i kk]))+100;

        a = ahx(1, ii);

%             scatter(a, X(:, i), X(:, kk), [], IDs, 'x', 'linewidth', 2);
        P1 = [XLIMS{k, i}(1) YLIMS{k, i}(1)];
        P2 = -P1;
        rectangle('Position', [P1 P2], 'FaceColor', [.9 .6 .6], 'LineStyle', 'none', 'parent', a)  %[x,y,w,h]
        scatter(a, X(idx3, i), X(idx3, kk), [], [.7 .7 .7], '.')
        scatter(a, X(idx1, i), X(idx1, kk), [], [.3 .3 .9], '.')
        scatter(a, X(idx2, i), X(idx2, kk), [], 'k', '.')

        plot(a, [mii maa], [mii maa],  '--c', 'linewidth', 2)
        plot(a, [mii maa] + BotmConsts([i i]), [mii maa] + BotmConsts([kk kk]),  '--m', 'linewidth', 2)
        if k==2
            xlabel(a, ['BOTM Template ' num2str(i)]);
        else
%                 set(a, 'xticklabel', []);
        end
        if kk==1 && i==1
            mysort.plot.legend(a, {{'.', 'color', [.3 .3 .9]}, ...
                                   {'.k'}, ...
                                   {'.', 'color', [.7 .7 .7]},...
                                   {'-', 'linewidth', 8, 'Color', [.9 .6 .6]}
                                   }, ...
                     {'Target X-Axis', 'Target Y-Axis', 'Other', 'Noise Region'}, 'SouthEast')
        end
        ylabel(a, ['BOTM Template ' num2str(kk)]);
        set(a, 'xlim', XLIMS{k, i}, 'ylim', YLIMS{k, i});
    end
end
% mysort.plot.figureTitle('');
% mysort.plot.figureDate();
mysort.plot.savefig(gcf, 'FigureForReviewer2_B', 'fig', 0, 'ai', 1);
mysort.plot.savefig(gcf, 'C:\Users\frankef\Dropbox\PaperBOTM1\SubmissionJCompNeuro\Revision1\FigureForReviewer2_B', 'fig', 0, 'ai', 1);

return

%%
XX = {TMbotm TMmatc};

BotmConsts = TMbotm(1,:) - TMmatc(1,:);

XLIMS = {[-90 110] [-10 80] [-150 200]
         [-150 200] [-30 100] [-150 200]};
YLIMS = {[-30 90] [-30 100] [-110 120]
         [-150 200] [-120 180] [-30 100]};
for k=1:1 %length(XX)
    X = XX{k};
    MI = min(X);
    MA = max(X);

    mysort.plot.figure([1000 550]);
    nR = 2; nC = 3; 
    ahx = mysort.plot.subplots([nR nC], 'matrix', 1, 'spacerX', .1);
    p = 1;
    for i=1:3
        idx1 = IDs==i; %ismember(IDs, [1 2]);
        k = 1;
        for kk=1:3
            if kk==i
                continue
            end
            idx2 = IDs == kk;
            idx3 = ~idx1 & ~idx2;
            
            mii= min(MI([i kk]));
            maa= max(MA([i kk]))+100;
            
            a = ahx(k, i);
            
%             scatter(a, X(:, i), X(:, kk), [], IDs, 'x', 'linewidth', 2);
            P1 = [XLIMS{k, i}(1) YLIMS{k, i}(1)];
            P2 = -P1;
            rectangle('Position', [P1 P2], 'FaceColor', [.9 .6 .6], 'LineStyle', 'none', 'parent', a)  %[x,y,w,h]
            scatter(a, X(idx3, i), X(idx3, kk), [], [.7 .7 .7], '.')
            scatter(a, X(idx1, i), X(idx1, kk), [], [.3 .3 .9], '.')
            scatter(a, X(idx2, i), X(idx2, kk), [], 'k', '.')

            plot(a, [mii maa], [mii maa],  '--c', 'linewidth', 2)
            plot(a, [mii maa] + BotmConsts([i i]), [mii maa] + BotmConsts([kk kk]),  '--m', 'linewidth', 2)
            if k==2
                xlabel(a, ['BOTM Template ' num2str(i)]);
            else
%                 set(a, 'xticklabel', []);
            end
            if k==1 && i==1
                mysort.plot.legend(a, {{'.', 'color', [.3 .3 .9]}, ...
                                       {'.k'}, ...
                                       {'.', 'color', [.7 .7 .7]},...
                                       {'-', 'linewidth', 8, 'Color', [.9 .6 .6]}
                                       }, ...
                         {'Target X-Axis', 'Target Y-Axis', 'Other', 'Noise Region'}, 'SouthEast')
            end
            ylabel(a, ['BOTM Template ' num2str(kk)]);
            set(a, 'xlim', XLIMS{k, i}, 'ylim', YLIMS{k, i});
            
%             set(a, 'xscale', 'log', 'yscale', 'log')
            k = k+1;
        end
    end
end
% mysort.plot.figureTitle('');
% mysort.plot.figureDate();
mysort.plot.savefig(gcf, 'FigureForReviewer2', 'fig', 0, 'ai', 1);
mysort.plot.savefig(gcf, 'C:\Users\frankef\Dropbox\PaperBOTM1\SubmissionJCompNeuro\Revision1\FigureForReviewer2', 'fig', 0, 'ai', 1);
