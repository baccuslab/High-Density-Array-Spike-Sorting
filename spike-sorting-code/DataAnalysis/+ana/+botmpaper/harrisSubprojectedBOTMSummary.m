    %% INIT
    close all

[REJECTIONPERF, CLASSPERF, DETPERF, CPUTIMEPERF, correctSpikesIdx, otherSpikesIdx, L, nTemplatesPerDataSet] = ana.botmpaper.loadHarrisEval();
%%
    mysort.plot.figure('w',1100, 'h', 800);
    ah = subplot(5,1,1); hold on
    mDP = squeeze(mean(DETPERF));
    sDP = squeeze(std(DETPERF));
    handles = matlabfilecentral.barweb.barweb(mDP, sDP, [], [], [], [], [], bone, [], L.methods(:,4), 1, 'axis');
    ylabel(ah(1), 'Det Perf [%]', 'fontsize', 12);
    set(ah(1), 'ylim', [0 100]);    
    
    ah(2) = subplot(5,1,2); hold on
    mCP = squeeze(mean(CLASSPERF));
    sCP = squeeze(std(CLASSPERF));    
    handles = matlabfilecentral.barweb.barweb(mCP, sCP, [], [], [], [], [], bone, [], [], 1, 'plot');
    set(ah(2), 'ylim', [0 101]);
    ylabel(ah(2), 'Class Perfo [%]', 'fontsize', 12);

    ah(3) = subplot(5,1,3); hold on
    mRP = squeeze(mean(REJECTIONPERF));
    sRP = squeeze(std(REJECTIONPERF));
    handles = matlabfilecentral.barweb.barweb(mRP, sRP, [], [], [], [], [], bone, [], [], 1, 'plot');
%     bar(ah(3), mRP, 'linewidth', 2)
    ylabel(ah(3), 'Rej Perf [%]', 'fontsize', 12);
    set(ah(3), 'ylim', [min(mRP(:)) 100]);
    set(ah(3), 'xticklabel', L.datasetNames(correctSpikesIdx))
    
     ah(4) = subplot(5,1,4); hold on
     PP = (REJECTIONPERF + CLASSPERF + DETPERF)/3;
    mRP = squeeze(mean(PP));
    sRP = squeeze(std(PP));
    handles = matlabfilecentral.barweb.barweb(mRP, sRP, [], [], [], [], [], bone, [], [], 1, 'plot');
%     bar(ah(3), mRP, 'linewidth', 2)
    ylabel(ah(4), 'Total Perf [%]', 'fontsize', 12);
    set(ah(4), 'ylim', [min(mRP(:)) 100]);
    set(ah(4), 'xticklabel', L.datasetNames(correctSpikesIdx))   
    
     ah(5) = subplot(5,1,5); hold on
     PP = CPUTIMEPERF;
    mRP = squeeze(mean(PP));
    sRP = squeeze(std(PP));
    handles = matlabfilecentral.barweb.barweb(mRP, sRP, [], [], [], [], [], bone, [], [], 1, 'plot');
    ylabel(ah(5), 'CPU Time [s]', 'fontsize', 12);
    set(ah(5), 'xticklabel', L.datasetNames(correctSpikesIdx))      
%     mysort.plot.savefig(gcf, 'BigPerformanceAllProcessingPlotHarris', 'fig', 0)
    
%     legend(ah(1), L.methods(:,4), 'location', 'northeastoutside')
%     legend(ah(2), L.methods(:,4), 'location', 'northeastoutside')
%     legend(ah(3), L.methods(:,4), 'location', 'northeastoutside')

%         tot = .5*squeeze(mean(DETPERF)) + .5*squeeze(mean(CLASSPERF));
    %% SAME PLOT AS ABOVE BUT GROUPED BY METHOD, NOT DATA
    preprocessingnames = L.datasetNames(otherSpikesIdx);
    interesing_preprocessing_idx = [1 4 5];
    int_prep_names = preprocessingnames(interesing_preprocessing_idx);
    int_prep_names
    % DETPERF is data sets x preprocessing x methods 
    [nDS nPreproc nMethods] = size(DETPERF);
    mysort.plot.figure('w',1100, 'h', 800);
    nR = 5; nC = nMethods; nP = nC*nR; p=0; ah = [];
    for m=1:nC
        cidx = 0;
        p=p+1; ah(p) = subplot(nR, nC, (nC*cidx)+m); hold on
        P = DETPERF(:,interesing_preprocessing_idx,m);
        mP = squeeze(mean(P));
        sP = squeeze(std(P));
        handles = matlabfilecentral.barweb.barweb(mP', sP', [], [], [], [], [], bone, [], [], 1, 'plot');
        title(L.methods(m,4));
        if m==1
            ylabel(ah(p), 'Det Perf [%]', 'fontsize', 12);
            set(ah(p), 'ylim', [0 0.4]);
        elseif m==2
            set(ah(p), 'ylim', [0 0.4]);
            set(ah(p), 'yticklabel', []);
        elseif m==3
            set(ah(p), 'ylim', [99.9 100.03]);
        else
            set(ah(p), 'ylim', [99.9 100.03]);
            set(ah(p), 'yticklabel', []);
        end
        
        cidx = cidx +1;
        p=p+1; ah(p) = subplot(nR, nC, (nC*cidx)+m); hold on
        P = CLASSPERF(:,interesing_preprocessing_idx,m);
        mP = squeeze(mean(P));
        sP = squeeze(std(P));
        handles = matlabfilecentral.barweb.barweb(mP', sP', [], [], [], [], [], bone, [], [], 1, 'plot');
        set(ah(p), 'ylim', [0 105]);
        if m==1
            ylabel(ah(p), 'Class Perf [%]', 'fontsize', 12);
        else
            set(ah(p), 'yticklabel', []);
        end 
        
        cidx = cidx +1;
        p=p+1; ah(p) = subplot(nR, nC, (nC*cidx)+m); hold on
        P = REJECTIONPERF(:,interesing_preprocessing_idx,m);
        mP = squeeze(mean(P));
        sP = squeeze(std(P));
        handles = matlabfilecentral.barweb.barweb(mP', sP', [], [], [], [], [], bone, [], [], 1, 'plot');
        set(ah(p), 'ylim', [50 105]);
        if m==1
            ylabel(ah(p), 'Rejec Perf [%]', 'fontsize', 12);
        else
            set(ah(p), 'yticklabel', []);            
        end   
        
        cidx = cidx +1;
        p=p+1; ah(p) = subplot(nR, nC, (nC*cidx)+m); hold on
        P = (REJECTIONPERF(:,interesing_preprocessing_idx,m) + CLASSPERF(:,interesing_preprocessing_idx,m) + DETPERF(:,interesing_preprocessing_idx,m))/3;
        mP = squeeze(mean(P));
        sP = squeeze(std(P));
        handles = matlabfilecentral.barweb.barweb(mP', sP', [], [], [], [], [], bone, [], [], 1, 'plot');
        set(ah(p), 'ylim', [0 105]);
        if m==1
            ylabel(ah(p), 'Perf [%]', 'fontsize', 12);
        else
            set(ah(p), 'yticklabel', []);            
        end    
        
        cidx = cidx +1;
        p=p+1; ah(p) = subplot(nR, nC, (nC*cidx)+m); hold on
        P = CPUTIMEPERF(:,interesing_preprocessing_idx,m);
        mP = squeeze(mean(P));
        sP = squeeze(std(P));
        handles = matlabfilecentral.barweb.barweb(mP', sP', [], [], [], [], [], bone, [], [], 1, 'plot');
        set(ah(p), 'ylim', [0 0.3]);
        if m==1
            ylabel(ah(p), 'CPU [s]', 'fontsize', 12);
%             
%         elseif m==2
%             set(ah(p), 'ylim', [0 0.3]);
%             set(ah(p), 'yticklabel', []);
%         elseif m==3
%             set(ah(p), 'ylim', [0 .008]);
%         else
%             set(ah(p), 'ylim', [0 .008]);
%             set(ah(p), 'yticklabel', []);
        end        
    end
%     mysort.plot.savefig(gcf, 'BigPerformancePlotHarris', 'fig', 0)
    
    %% Make CPU Plot
    P = CPUTIMEPERF(:,interesing_preprocessing_idx,:);
    mP = squeeze(mean(P));
    sP = squeeze(std(P));    
    mysort.plot.figure('w',1100, 'h', 800);
    ah = axes('fontsize', 14);
    bar(mP');
    set(ah, 'xtick', 1:size(L.methods(:,4),1), 'xticklabel', L.methods(:,4));
    legend(int_prep_names, 'location', 'northeast', 'fontsize', 14);
    ylabel('CPU Time [s]', 'fontsize', 14);
%     set(ah, 'YScale', 'log');
%     mysort.plot.savefig(gcf, 'CPUPlotHarris', 'fig', 0)

    %% 
    mysort.plot.figure('w',1100, 'h', 800);
    ah = subplot(1,2,1);
    