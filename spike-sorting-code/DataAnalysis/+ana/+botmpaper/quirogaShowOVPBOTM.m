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
%%
b = 1;
n=1; 

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
        
        % Noise estimation
        spikeEpochs = mysort.epoch.merge([tgdf(:,2)-20 tgdf(:,2)+80]);
        noiseEpochs = mysort.epoch.flip(spikeEpochs, size(GT.X,2));
        noiseEpochs = mysort.epoch.removeShort(noiseEpochs, Tf);
        DS = mysort.ds.Matrix(GT.X', 32000, 'Q', 1, 1);
        noiseSnippets = DS.getWaveform(noiseEpochs(:,1), 0, Tf);
        nNoise = size(noiseSnippets,1);
        
        XCC = mysort.noise.XCorrContainer(DS, Tf-1, 'noiseEpochs', noiseEpochs);
        Cest = XCC.getCte4Channels(1:size(DS,2));        
        Cest = .5*Cest + .5*diag(diag(Cest));
         
        Covest = struct();
        Covest.CCol = mysort.noise.Cte2Ccol(Cest, size(DS,2));
        
        GTSortingContainer = mysort.spiketrain.SpikeSortingContainer('gt', tgdf, 'wfDataSource', DS);
        DS.addSpikeSorting(GTSortingContainer);
        
        
	tT = mysort.wf.v2t(T, 1);
	[tTA tau] = mysort.wf.tAlignOnCorrelation(tT, 'debug', 1, 'trunc', 1);
    vTA = mysort.wf.t2v(tTA);
    
    figure;
    subplot(2,1,1)
    plot(T')
    subplot(2,1,2)
    plot(vTA')
    
    botm = mysort.sorters.BOTM(Covest, tTA, 'adaptOnInit', 1);
    botm.DH = DS;
    botm.sortMatrix(DS(1:1000,:)');    
    botm.name = 'BOTM';
      
    botmOvp = mysort.sorters.BOTMOvp(Covest, tTA, 'adaptOnInit', 1);
    botmOvp.DH = DS;
    botmOvp.sortMatrix(DS(1:1000,:)');
    botmOvp.name = 'BOTMOvp';
    
    mysort.plot.SliderDataAxes({DS, botm, botmOvp}, 'channelSpacers', [0,  0, 0])

    
%     DS3 = mysort.ds.Matrix(Xfil', D.srate, 'bla');
%     BOTMSortingContainer1 = mysort.spiketrain.SpikeSortingContainer('gt', botmgdf1, 'wfDataSource', DS);
%     DS3.addSpikeSorting(BOTMSortingContainer1);
% 
%     DS4 = mysort.ds.Matrix(Xfil', D.srate, 'bla');
%     BOTMSortingContainer2 = mysort.spiketrain.SpikeSortingContainer('gt', botmgdf2, 'wfDataSource', DS);
%     DS4.addSpikeSorting(BOTMSortingContainer2);
%     
%     mysort.plot.SliderDataAxes({DS1, DS2, DS3, botm}, 'channelSpacers', [400, 400, 400, 0])