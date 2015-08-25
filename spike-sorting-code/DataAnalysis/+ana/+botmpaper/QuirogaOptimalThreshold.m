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
savefig__ = 0;



lambda = .5;
bUseCest = 0;
figSavePath = 'C:\Dropbox\PaperBOTM1\DetectorFigures\';
suffix = '';

% lambda = 0;
% bUseCest = 1;
% figSavePath = 'C:\Users\frankef\Dropbox\PaperBOTM1\DetectorFigures_Cest\';
% suffix = 'Cest_lamb0';

%%
RR = [];
for b = 1:4 %length(benchmarks)
    for n = 1:4 %length(noiselevels)
        savePath = fullfile(dpath, benchmarks{b}, noiselevels{n});
        if ~exist(savePath, 'file')
            mkdir(savePath);
        end      
    
        pic = 0;
        fprintf('Starting with %s %s\n', benchmarks{b}, noiselevels{n});
        close all
        quirogaFile = sprintf('C_%s_noise%s.mat', benchmarks{b}, noiselevels{n});        
        
        bufferFile = fullfile(savePath, ['botmOutputDistributionsThrehold' suffix '.mat']);
        if exist(bufferFile, 'file')
            load(bufferFile);
        else
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
            Tf = size(T,2);

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
            if bUseCest==1
                C_ = (1-lambda)*Cest + lambda*targetC;
            else
                C_ = (1-lambda)*C + lambda*targetC;
            end

            % Filtering of the data and estimation of expected peaks
            Y = zeros(size(GT.X));
            Yxc = [];
            for i=1:3
                t = T(i,:);
                f = t/C_;
                Y(i,:)  = conv(GT.X, flipud(f(:)), 'same') - .5*t*f';
                Yxc(i,:) = conv(t, flipud(f(:)), 'same');
            end

            if 0
                figure
                plot(Y')
                hold on
                plot(tgdf(:,2), zeros(size(tgdf,1),1), 'cx', 'markersize', 14, 'linewidth',3 )
            end

            %% Find peaks
            pks = {};
            locs = {};
            peaks = {};
            for i=1:3
                [pks, locs] = findpeaks(Y(i,:), 'MINPEAKDISTANCE', 15);
                peaks{i} = [locs(:) pks(:)];
            end

            %%
            save(bufferFile, 'peaks', 'C', 'C_', 'T', 'tgdf', 'cutLeft', 'Tf');
        end
        %% annotation of peaks
        OVERLAP_DIST = 15;
        OTHERSPIKE_DIST = 6;
        SINGLESPIKE_DIST = 3;
        D = [];
        for i=1:3
            p = peaks{i};
            p(end, 3) = 0;
            mySpikeTimes = tgdf(tgdf(:,1)==i,2) + cutLeft + 1;
            otherSpikeTimes = tgdf(tgdf(:,1)~=i,2) + cutLeft + 1;
            for k=1:size(p,1)
                dists = abs(mySpikeTimes - p(k,1));
                di = find(dists < SINGLESPIKE_DIST);
                do = abs(otherSpikeTimes - p(k,1));
                bIsOvp = any(do < OVERLAP_DIST);
                bIsOtherSpike = any(do < OTHERSPIKE_DIST);
                if ~isempty(di);
                    D(end+1) = dists(di);
                    if bIsOvp
                        p(k,3) = 2; % ovp
                    else
                        p(k,3) = 1; % single spike
                    end
                else
                    if bIsOtherSpike
                        p(k,3) = 3; % other spike
                    else
                        p(k,3) = 4; % noise
                    end
                end
            end
            peaks{i} = p;
        end
        
        %% Compute Error over threhold
        
%         pNRange = .5:.001:.999999;
%         thresholdRange = log(pNRange);
%         TP = zeros(3, length(thresholdRange));
%         FN = zeros(3, length(thresholdRange));
%         FP = zeros(3, length(thresholdRange));
%         CL = zeros(3, length(thresholdRange));
%         allPeaks = [];
%         for i=1:3
%             p = peaks{i};
%             allPeaks(:,i) = p(:,2);
%         end
%         [allMaxPeaks allPeakClassifications] = max(allPeaks,2);
%         for i=1:3
%             mySingleSpikes = p(:,3) == 1;
%             myOvpSpikes    = p(:,3) == 2;
%             myOtherSpikes  = p(:,3) == 3;
%             myNoise        = p(:,3) == 4;
%             
%             for t=1:length(thresholdRange)
%                 thr = thresholdRange(t);
%                 
%                 TP(i,t) = sum( p(mySingleSpikes,2)>thr );
%                 FP(i,t) = sum( p(myNoise,2)>thr);
%                 FP(i,t) = sum( p(myNoise,2)>thr);
%             end
%         end
        
        if 0
            %%
            GT = ana.botmpaper.E10_preprocessing(dpath, quirogaFile, initCutLeft, initTf, dpath);
            gt_spikes = mysort.epoch.extractWaveform(GT.X, [GT.gdf(:,2)-initCutLeft GT.gdf(:,2)-initCutLeft+initTf-1]);
            % Filtering of the data and estimation of expected peaks
            Y = zeros(size(GT.X));
            Yxc = [];
            for i=1:3
                t = T(i,:);
                f = t/C_;
                Y(i,:)  = conv(GT.X, flipud(f(:)), 'same') - .5*t*f';
                Yxc(i,:) = conv(t, flipud(f(:)), 'same');
            end
        end
        if 0
            %%            
            seli = 1;
            figure
            plot(Y(seli,:));
            hold on
            plot(tgdf(tgdf(:,1)==seli,2), zeros(size(tgdf(:,1)==seli,2),1), 'cx', 'markersize', 14, 'linewidth',3 )
            plot(tgdf(tgdf(:,1)~=seli,2), zeros(size(tgdf(:,1)~=seli,2),1), 'ro', 'markersize', 14, 'linewidth',3 )
        end
        %% 
        fh2 = mysort.plot.figure('w',1200, 'h', 500);
        ph=0; ah=[]; nR= 4; nC= 3;
        cols = {'g', 'b', 'm', 'r'};
        names = {'Single Spikes', 'Overlaps', 'Other Spikes', 'Noise'};
        thr1 = log(1-3*0.1);
        thr2 = log(1-3*0.0000001);
        % compute min max values for each neuron
        mima = [];
        for i=1:3
            p=peaks{i};
            mima(i,:) = [min(p(:,2)) max(p(:,2))];
        end
        for j=1:4
            for i=1:3
                p=peaks{i};
                
                x = p(round(p(:,3))==j,2);                
                centers = linspace(mima(i,1), mima(i,2), 100);
                edges = centers + (centers(2)-centers(1))/2;
        
                n1 = histc(x, edges);
                
                ph=ph+1; ah(ph) = subplot(nR, nC, ph);
                hold on
                bar(centers, n1(1:end), 'facecolor', cols{j}, 'EdgeColor', 'none', 'BarWidth', 1);
%                 if j==1
%                     title([names{j} ' - Neuron ' num2str(i)])
%                 else
%                     title(names{j})
%                 end
                set(ah(ph), 'xlim', [centers(1) centers(end)], 'ylim', [0 max(n1)]);
                plot(thr1*[1 1], [0 max(n1)], '--', 'color', [.7 .7 .7], 'linewidth', 1);
                plot(thr2*[1 1], [0 max(n1)], '--', 'color', [.5 .5 .9], 'linewidth', 1);
%                 legend(['N = ' num2str(length(x))])
                if j==4
                    xlabel('local maxima in BOTM Output');
                else
                    set(gca, 'xticklabel', []);
                end
                if i==1
                    ylabel('Count');
                end
            end
        end
        mysort.plot.figureTitle(sprintf('%s %s', benchmarks{b}, noiselevels{n}));
        filePrefix = sprintf('%s_%s', benchmarks{b}, noiselevels{n});
        mysort.plot.savefig(fh2, fullfile(figSavePath, [filePrefix 'detector_distributions_threshold']), 'ai', 1, 'fig', 0);

        fprintf('Done.\n');
    end
end
