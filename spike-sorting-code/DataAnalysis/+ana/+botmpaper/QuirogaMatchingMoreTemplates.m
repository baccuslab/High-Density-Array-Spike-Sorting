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
RRRR = {};
for b = 1:4 %length(benchmarks)
    for n = 1:4 %length(noiselevels)
        savePath = fullfile(dpath, 'MoreTemplateAnalysis', benchmarks{b}, noiselevels{n});
        if ~exist(savePath, 'file')
            mkdir(savePath);
        end      
        bufferFile = fullfile(savePath, ['QuirogaMatchingMoreTemplatesResults_b' num2str(b) '_n' num2str(n) '.mat']);
        if ~exist(bufferFile, 'file')
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
            targetCest = dCest;        

            lambda = .5;
            Cest_ = (1-lambda)*Cest + lambda*targetCest; 
            T_org = T;

            parfor ti = 1:20
                for rep = 1:5
                    close all

                    %%
                    T = ana.botmpaper.QuirogaMatchingMoreTemplatesHelper(T_org, ti);
                    figure
                    plot(T(1:3,:)', 'linewidth', 2);
                    hold on
                    plot(T(4:end,:)', 'k', 'linewidth', 1);
                    %%


                    DATASETS = {gt_spikes, T, Cest, 'corr, Cest', 0 , [] 
                                noiseSnippets,  T, Cest, 'N, Cest', 0 , [] 
                                gt_spikes, T, dCest, 'corr, dC', 0 , [] 
                                noiseSnippets,  T, dCest, 'N, dC', 0 , []   
                                gt_spikes, T, Cest_, 'corr, Cest_', 0 , [] 
                                noiseSnippets,  T, Cest_, 'N, Cest_', 0 , [] 
                                };

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

                    noiseCorrectPairs = repmat([1  2], size(DATASETS,1)/2, 1) + repmat((0:(size(DATASETS,1)/2 -1))'*2, 1, 2);

                    %% Use Classification responses to noise and to the correct template to compute the
                    % detection performance
                    detThr = zeros(size(noiseCorrectPairs,1), nMethods);
                    detPerf = detThr; 

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

                        end
                    end  
                    correctSpikesIdx = 1:2:size(DATASETS,1);

    %                 %% Make big performance plot, all vs. all
    %                 mysort.plot.figure('w',1000,'h',1000);
    %                 p = 1;
    %                 ah = [];
    %                 
    %                 for d=correctSpikesIdx
    %                     lT = DATASETS{d,2};
    %                     nT  = size(lT,1); 
    %                     nS  = size(DATASETS{d,1},1);
    %                     for m=1:nMethods
    %                         ah(p) = subplot(length(correctSpikesIdx),nMethods,p);
    %                         bar(1:nT, RES(d).classifications(1:nT,m))
    %                         title(methods{m,4});
    %                         if m == 1
    %                             ylabel(DATASETS{d, 4})
    %                         end
    %                         p=p+1;
    %                     end
    %                 end
    %                 set(ah, 'ylim', [0 nS+10]);
    %                 pic=pic+1; mysort.plot.savefig(gcf, fullfile(savePath, [sprintf('%03d_', pic) 'classifications']));        

                    % Make big performance plot, only performance
                    ah = [];
                    PERF = [];
                    for di=1:length(correctSpikesIdx)
                        d = correctSpikesIdx(di);
                        for m=1:nMethods
                            PERF(di,m) = 100*RES(d).performance(m);
                        end
                    end
                    tot = (-100*detPerf + PERF )/2;                

    %                 %
    %                 mysort.plot.figure('w',1100, 'h', 800);
    %                 ah = subplot(3,1,1);
    %                 ah(2) = subplot(3,1,2);
    %                 ah(3) = subplot(3,1,3);
    % 
    %                 set(ah, 'fontsize', 12);
    %                 bar(ah(1), -100*detPerf)
    %                 ylabel(ah(1), 'Detection Perf [%]', 'fontsize', 12);
    %                 set(ah(1), 'ylim', [60 100]);    
    % 
    %                 bar(ah(2), PERF);    
    %                 set(ah(2), 'ylim', [0 100]);
    %                 ylabel(ah(2), 'Class Perf [%]', 'fontsize', 12);
    % 
    % 

    %                 bar(ah(3), tot);    
    %                 set(ah(3), 'ylim', [0 100]);
    %                 ylabel(ah(3), 'Total [%]', 'fontsize', 12);  
    %                 set(ah(3), 'ylim', [min(tot(:)) 100]);
    %                 set(ah(3), 'xticklabel', DATASETS(correctSpikesIdx,4))
    % 
    %                 legend(ah(1), methods(:,4), 'location', 'northeastoutside')
    %                 legend(ah(2), methods(:,4), 'location', 'northeastoutside')    
    %                 legend(ah(3), methods(:,4), 'location', 'northeastoutside')
    %                 mysort.plot.figureTitle([benchmarks{b} ' ' noiselevels{n}])
    %                 pic=pic+1; mysort.plot.savefig(gcf, fullfile(savePath, [sprintf('%03d_', pic) 'final_results']));        

                    %%
                    R = struct();
                    R.detPerf = detPerf;
    %                 R.RES = RES;
                    R.PERF = PERF;
                    R.tot = tot;
                    R.T = T;
                    RRR{ti, rep} = R;
            %         clear GT X DATASETS noiseSnippets XCC meanFreeNoise SP NP mx mxidx K KS KN otherSpikes
            % %         save(fullfile(savePath, 'final_results'))
            %         save(fullfile(savePath, 'final_results_small'), 'tgdf', 'RES', 'detPerf', 'methods', 'datasetNames');
                    fprintf('Done with %s\n', [benchmarks{b} ' ' noiselevels{n}]);
                end
            end
            save(bufferFile, 'RRR', 'methods', 'datasetNames');
        else
            load(bufferFile);
        end
        RRRR{b, n} = RRR;
    end
end
save('QuirogaMatchingMoreTemplatesResults.mat', 'RRRR', 'methods', 'datasetNames');

%%
Cest_Idx = 3;
botmIdx = 5;
X = zeros(4, 4, 20, 5);
for b=1:4
    for n=1:4
        for ti=1:20
            for rep = 1:5
                X(b,n,ti,rep) = RRRR{b, n}{ti,rep}.tot(Cest_Idx,botmIdx);
            end
        end
    end
end
XM_b = squeeze(mean(mean(X,4),1))';
XS_b = squeeze(std(std(X,[],4),1))';
XM_n = squeeze(mean(mean(X,4),2))';
XS_n = squeeze(std(std(X,[],4),1))';
save('QuiroagaMatchingMoreTemplatesCondensedResults.mat', 'XM_b', 'XS_b', 'XM_n', 'XS_n');
%%
b = 4;
n = 4;
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
%%
load('QuiroagaMatchingMoreTemplatesCondensedResults.mat', 'XM_b', 'XS_b', 'XM_n', 'XS_n');

T1 = ana.botmpaper.QuirogaMatchingMoreTemplatesHelper(T, 16);

mysort.plot.figure([1000 600])
subplot(2,5,1:3)
% plot(XM_b)
boundedline(1:20, XM_b, cat(2, permute(XS_b,[1 3 2]), permute(XS_b,[1 3 2])))
title('Mean Over Benchmarks')
hl = legend('Easy1', 'Easy2', 'Difficult1', 'Difficult2'); set(hl, 'box', 'off');
ylabel('Total Performance [%%]')
        
subplot(2,5,6:8)
boundedline(1:20, XM_n, cat(2, permute(XS_n,[1 3 2]), permute(XS_n,[1 3 2])))
title('Mean Over Noise')
hl = legend('N .05', 'N .10', 'N .15', 'N .2'); set(hl, 'box', 'off');
xlabel('Additional Templates');
ylabel('Total Performance [%%]')



subplot(2,5, 4:5); hold on
plot(T1(1:3,:)', 'linewidth', 2);
axis tight
box off
ylabel('Amplitude')
% xlabel('Time [samples]')
title('Original Templates Difficult2');

subplot(2,5, 9:10); hold on
plot(T1(1:3,:)', 'm', 'linewidth', 2);
plot(T1(4:end,:)', 'c', 'linewidth', 1);
axis tight
box off
ylabel('Amplitude')
xlabel('Time [samples]')
title('Manipulated Template set (+16)');
% mysort.plot.savefig(gcf, 'QuirogaMatchingMoreTemplates', 'ai', 1, 'fig', 0);
% mysort.plot.savefig(gcf, 'C:\Users\frankef\Dropbox\PaperBOTM1\SubmissionJCompNeuro\Revision1\QuirogaMatchingMoreTemplates', 'ai', 1, 'fig', 0);