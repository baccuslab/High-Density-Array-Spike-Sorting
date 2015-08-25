function handles = PlotPSTHButton(hObject, eventdata, handles)
    % hObject    handle to PlotAmplitudeDriftButton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    [P PP] = dataviewer.util.getHandleValues(handles);  
    EsTl = handles.DH.getEventSampleAndTrialLength('trialIDs', P.trialIDs, 'eventIDs', P.eventIDs);
    if isempty(EsTl)
        handles.warning_func('Not trials with matching constraints!');
    end
    validTrialIDs = unique(EsTl(:,1));
    invalidTrials = setdiff(P.trialIDs, validTrialIDs);
    GDF  = handles.DH.getSpikeTrain('otherP', P, 'trialIDs', validTrialIDs, 'alignOnEvent', P.eventIDs);    
    if isempty(EsTl)
        handles.warning_func('Not trials with matching constraints!');
    end
    srate = handles.DH.getSamplesPerSecond();
    
    trialIDs = EsTl(:,1);
    trialIdx = EsTl(:,2);
    trialLen = EsTl(:,3);
    eventSam = EsTl(:,4);  
    
    nValidTrials = length(trialIDs);
    
    if isempty(P.from)
        P.from = 1;
    end
    if isempty(P.to)
        P.to = max(trialLen);
    end

    % Event Stats
    minEvSample  = min(eventSam);
    maxEvSample  = max(eventSam);
    medianEvSample = median(eventSam);
    stdEvSample  = sqrt(eventSam);
    
    % sample Stats
    minSample = min(GDF(:,5));
    maxSample = max(GDF(:,5));
    
    % Binnings  
    binSizeInMS = P.binsize;
    binSize     = binSizeInMS * srate/1000;
    scalingFactor = srate/binSize;
    
    edges = minSample:binSize:maxSample;
    % ignore first bin
    edges = edges(2:end);
    nBins = length(edges)-1;
    binCenters = (edges(1:end-1)+edges(2:end))/2;
    if nBins == 0
        handles.warning_func('No valid bins to calculate PSTH in!');
        return
    end
    
    % Reference Bins
    referenceSamples = [minSample minSample+srate+1000];   % For BG Rate
%     if P.from >= referenceSamples(2) || P.to <= referenceSamples(1)
%         referenceBins = [];
%     else
        referenceBins = find(edges<referenceSamples(2) & edges>referenceSamples(1));
        referenceBins = referenceBins(1:end-1);
%     end

    % build big binning matrix
    COUNTS = zeros(nValidTrials, nBins);
    for i=1:nValidTrials
        tmp = histc(GDF(GDF(:,1)==trialIDs(i),5), edges);
        COUNTS(i,:) = tmp(1:end-1);  % the last value of histc counts equality with edges(end)
    end

    % estimate the number of valid trials per bin
    validTrialsPerBin = zeros(1, nBins);
    variancePerBin    = zeros(1, nBins);
    standardErrorOfTheMean = zeros(1, nBins);
    for i=1:nBins
        validTrialsPerBin(i) = sum(trialLen-eventSam>edges(i+1));
        variancePerBin(i)    = var(COUNTS(trialLen-eventSam>edges(i+1), i));
        standardErrorOfTheMean(i) = sqrt(variancePerBin(i)) / sqrt(validTrialsPerBin(i));
    end    
    
    % standard error of the mean
    counts = sum(COUNTS,1);
    meanCounts = counts./validTrialsPerBin;
    meanRate = meanCounts * scalingFactor;
    handles.new_figure_handles = mysort.plot.figure();
    plot(binCenters/srate, meanRate , 'k','Linewidth',2);
    hold on
    % debug:
    % figure; plot(standardErrorOfTheMean);hold on;plot(variancePerBin,'g');plot(validTrialsPerBin,'r');legend('SEM', 'VAR', 'N');
    if length(referenceBins)>1       
        bgRate = mean(meanRate(referenceBins));
        plot(binCenters([1 end])/srate, bgRate*ones(1,2), ':g', 'linewidth', 2);
        plot(binCenters/srate, meanRate + 1.96*standardErrorOfTheMean* scalingFactor, '-', 'linewidth', 1, 'color', .8*ones(1,3));
        plot(binCenters/srate, meanRate - 1.96*standardErrorOfTheMean* scalingFactor, '-', 'linewidth', 1, 'color', .8*ones(1,3));
    end
    
    ylims = get(gca, 'ylim');
    set(gca, 'ylim', [0 ylims(2)]);
    plot([0 0],[0 ylims(2)],'Color', 'r','LineStyle',':','Linewidth',2);
    set(gca, 'xlim', [edges(2) edges(end)]/srate);
    xlabel(['tim in s (binsize: ' num2str(binSizeInMS) 'ms)'])
    ylabel('rate [Hz]');
    legend('E[rate]', 'BG rate', '+1.96*SEM', '-1.96*SEM', 'Stimulus');
    mysort.plot.figureTitle(['PSTH for ' PP.event.val{1} ' N = ' num2str(nValidTrials) ' (' P.figureTitle ')']);
    mysort.plot.figureName('PSTH'); 
   	handles.new_figure_handles = gcf;
    