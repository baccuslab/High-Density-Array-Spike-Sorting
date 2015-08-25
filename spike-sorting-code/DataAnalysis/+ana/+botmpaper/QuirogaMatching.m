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
for b = 1:4 %length(benchmarks)
    for n=1:4 %length(noiselevels)
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
        
        % Make plot of count histograms for noise and signal
        noiseCorrectPairs = repmat([1  2], size(DATASETS,1)/2, 1) + repmat((0:(size(DATASETS,1)/2 -1))'*2, 1, 2);
        mysort.plot.figure('w',1000,'h',1000); p = 0;
        correctTemplateIdx = 1;
        for pairi = 1:size(noiseCorrectPairs,1)
            correctDSIdx = noiseCorrectPairs(pairi,1);
            noiseDSIdx = noiseCorrectPairs(pairi,2);

%             neuronFR = 90; %Hz
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
        pic=pic+1; mysort.plot.savefig(gcf, fullfile(savePath, [sprintf('%03d_', pic) 'counts_noise_signal']));      
        
        %% Use Classification responses to noise and to the correct template to compute the
        % detection performance
        detThr = zeros(size(noiseCorrectPairs,1), nMethods);
        detPerf = detThr;
        fh1 = mysort.plot.figure('w',1000,'h',1000); p=0;
        fh2 = mysort.plot.figure('w',1200, 'h', 900);
        ah2    = subplot(3,2,3); hold on
        ah2(2) = subplot(3,2,4); hold on
        ah2(3) = subplot(3,2,5); hold on
        ah2(4) = subplot(3,2,6); hold on
        ah2(5) = subplot(3,2,1); hold on
        ah2(6) = subplot(3,2,2); hold on
        set(ah2, 'fontsize', 14);        

        for pairi = 1:size(noiseCorrectPairs,1)
            correctDSIdx = noiseCorrectPairs(pairi,1);
            noiseDSIdx = noiseCorrectPairs(pairi,2);

            neuronFR = 90; %Hz
            noisePriorWeight = samplesPerSecond/neuronFR;
%             noisePriorWeight = 1;
            nCorrect = size(RES(correctDSIdx).TM,1);
            nNoise   = size(RES(noiseDSIdx).TM,1);

            for m=1:nMethods
                noiseResponse = max(squeeze(RES(noiseDSIdx).TM(:,m,:)),[],2);
                signalResponse = max(squeeze(RES(correctDSIdx).TM(:,m,:)),[],2);
                noisefun = @(thr) sum(noiseResponse < thr)/nNoise;
                signafun = @(thr) sum(signalResponse < thr)/nCorrect;

                if strcmp(methods{m,3}, 'minDetection')
                    E = @(thr)  (noisePriorWeight*noisefun(thr) - signafun(thr))/noisePriorWeight;
                    thr0 = .8 * mean(noiseResponse);
                else
                    E = @(thr) (-noisePriorWeight*noisefun(thr) + signafun(thr))/noisePriorWeight;        
                    thr0 = .8 * mean(signalResponse);
                end
                if ~isempty(strfind(methods{m, 4}, 'Opt'))
                    detThr(pairi,m) = 0; 
                    detPerf(pairi, m) = E(0);
                else
                    [detThr(pairi,m) detPerf(pairi, m)] = fminsearch(E, thr0);
                end
                if ismember(m, [2 4]) && pairi == 7
                    c = 1;
                    if m==2
                        c = 0;
                    end
                    RRR = [];
                    mima = sort([-1*detThr(pairi,m) 3*detThr(pairi,m)]);
                    mima = [min(-50, mima(1)) max(100, mima(2))];
                    xr = linspace(mima(1),mima(2), 1000); %0:.05:10;                    
                    [nn,nxout] = hist(noiseResponse, 100);
                    [sn,sxout] = hist(signalResponse, 100);
%                     xr = -100:1:1000;
                    for kk=1:length(xr)
                        RRR(kk,1) = 1+ E(xr(kk));
                        RRR(kk,2) = noisefun(xr(kk));
                        RRR(kk,3) = signafun(xr(kk));
                    end
                    plot(ah2(5+c), nxout, nn*noisePriorWeight, 'r', 'linewidth', 3)
                    plot(ah2(5+c), sxout, sn, 'g', 'linewidth', 3)
                    ylabel(ah2(5+c), 'count', 'fontsize', 14);
                    title(ah2(5+c), methods{m, 4}, 'fontsize', 14);
                    
                    
                    plot(ah2(1+c), xr, RRR(:,2), 'r', 'linewidth', 3)
                    plot(ah2(1+c), xr, RRR(:,3), 'g', 'linewidth', 3)
                    ylabel(ah2(1+c), 'cumulative density', 'fontsize', 14);
                    set(ah2(1+c), 'ylim', [-.05 1.02]);
                    
                    plot(ah2(3+c), xr, RRR(:,1), 'b', 'linewidth', 3)
                    hold on
                    plot(ah2(3+c), detThr(pairi,m), 1+detPerf(pairi, m), 'gd', 'linewidth', 2);
                    ylabel(ah2(3+c), 'Error Function', 'fontsize', 14);
                    xlabel(ah2(3+c), 'Detector Treshold', 'fontsize', 14);
                    
                    if c == 1
                        h = legend(ah2(1+c), 'noise', 'signal');
                        set(h, 'location', 'southeast', 'fontsize', 14);
                        h = legend(ah2(3+c), 'Error Function', 'Opt. Threshold'); 
                        set(h, 'location', 'northeast', 'fontsize', 14);

                        linkaxes(ah2([1 3 5]), 'x');
                        linkaxes(ah2([2 4 6]), 'x');
                    end
                end                
                if 1
                    figure(fh1);
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
        pic=pic+1; mysort.plot.savefig(fh1, fullfile(savePath, [sprintf('%03d_', pic) 'detection_performance_estimation']));        
        pic=pic+1; mysort.plot.savefig(fh2, fullfile(savePath, [sprintf('%03d_', pic) 'detection_performance_for_maha_match']));            
        %% Make big performance plot, all vs. all
        mysort.plot.figure('w',1000,'h',1000);
        p = 1;
        ah = [];
        correctSpikesIdx = 1:2:size(DATASETS,1);
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
        pic=pic+1; mysort.plot.savefig(gcf, fullfile(savePath, [sprintf('%03d_', pic) 'classifications']));        
        
        % Make big performance plot, only performance
        ah = [];
        PERF = [];
        for di=1:length(correctSpikesIdx)
            d = correctSpikesIdx(di);
            for m=1:nMethods
                PERF(di,m) = 100*RES(d).performance(m);
            end
        end

        %
        mysort.plot.figure('w',1100, 'h', 800);
        ah = subplot(3,1,1);
        ah(2) = subplot(3,1,2);
        ah(3) = subplot(3,1,3);
 
        set(ah, 'fontsize', 12);
        bar(ah(1), -100*detPerf)
        ylabel(ah(1), 'Detection Perf [%]', 'fontsize', 12);
        set(ah(1), 'ylim', [60 100]);    

        bar(ah(2), PERF);    
        set(ah(2), 'ylim', [0 100]);
        ylabel(ah(2), 'Class Perf [%]', 'fontsize', 12);


        tot = (-100*detPerf + PERF )/2;
        bar(ah(3), tot);    
        set(ah(3), 'ylim', [0 100]);
        ylabel(ah(3), 'Total [%]', 'fontsize', 12);  
        set(ah(3), 'ylim', [min(tot(:)) 100]);
        set(ah(3), 'xticklabel', DATASETS(correctSpikesIdx,4))
        
        legend(ah(1), methods(:,4), 'location', 'northeastoutside')
        legend(ah(2), methods(:,4), 'location', 'northeastoutside')    
        legend(ah(3), methods(:,4), 'location', 'northeastoutside')
        mysort.plot.figureTitle([benchmarks{b} ' ' noiselevels{n}])
        pic=pic+1; mysort.plot.savefig(gcf, fullfile(savePath, [sprintf('%03d_', pic) 'final_results']));        
        %
        clear GT X DATASETS noiseSnippets XCC meanFreeNoise SP NP mx mxidx K KS KN otherSpikes
        save(fullfile(savePath, 'final_results'))
        save(fullfile(savePath, 'final_results_small'), 'tgdf', 'RES', 'detPerf', 'methods', 'datasetNames');
        fprintf('Done with %s\n', [benchmarks{b} ' ' noiselevels{n}]);
    end
end

        
  
    
% 
%         %% Make plot of theoretical and actual distribution of norms
%         mysort.plot.figure('w', 1200, 'h', 1000);
%         p = 1;
%         nP = size(noiseCorrectPairs,1);
%         ah = [];    
%         for pairi = 1:nP
%             correctDSIdx = noiseCorrectPairs(pairi,1);
%             noiseDSIdx = noiseCorrectPairs(pairi,2);
%             S      = DATASETS{correctDSIdx,1};
%             lT     = DATASETS{correctDSIdx,2};
%             Cpair  = DATASETS{correctDSIdx,3};
%             U = chol(Cpair);
%             N      = DATASETS{noiseDSIdx,1};
%     %         Cnoise = DATASETS{noiseDSIdx,3};
% 
%             dims = size(S,2);
%             corrT = repmat(mean(S), size(S,1), 1);
%             R = S - corrT;
%             R = R/U;
%             e = norm(corrT/U);
%             N = N/U;
% 
%             correct_distr = sum(R.^2,2);
%     %         correct_edges = linspace(min(correct_distr), max(correct_distr), 1000);
%             if pairi == nP
%                 correct_edges = linspace(0, 110, 100);
%             else
%                 correct_edges = linspace(0, 400, 100);
%             end
%             wcorr = correct_edges(2) - correct_edges(1);
%             correct_bincenters = correct_edges + wcorr/2;
%             correct_distr_counts = histc(correct_distr, correct_edges);
%             correct_distr_pdf    = correct_distr_counts/(sum(correct_distr_counts)*wcorr);
%             theo_correct = chi2pdf((correct_edges), dims);
% 
%             noise_distr = sum(N.^2,2);
%     %         noise_edges = linspace(min(noise_distr), max(noise_distr), 1000);
%             noise_edges = correct_edges;
%             wnoise = noise_edges(2) - noise_edges(1);
%             noise_bincenters = noise_edges + wnoise/2;        
%             noise_distr_counts = histc(noise_distr, noise_edges);
%             noise_distr_pdf    = noise_distr_counts/(sum(noise_distr_counts)*wnoise);
%             theo_noise  = theo_correct; %ncx2pdf(noise_edges, dims, e^2);
% 
%             ah(p) = subplot(2,nP,p); hold on
%             plot(ah(p), correct_bincenters, correct_distr_pdf, 'g', 'linewidth', 2);
%             plot(ah(p), correct_bincenters, theo_correct, 'k', 'linewidth', 2);
%             title(DATASETS{correctDSIdx,4});
% 
%             ah(p) = subplot(2,nP,p+nP); hold on
%             plot(ah(p), noise_bincenters, noise_distr_pdf, 'r', 'linewidth', 2);
%             plot(ah(p), noise_bincenters, theo_noise, 'k', 'linewidth', 2);
%             axis tight
%             p = p+1;
%         end
%     %     pic=pic+1; mysort.plot.savefig(gcf, fullfile(savePath, [sprintf('%03d_', pic) 'distributionComparisons']));

