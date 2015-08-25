function ah = xcorr_raw(spiketrains_or_gdf, varargin)
    P.IDs = [];
    P.T = [];
    P.maxLag = 50;    % in milliseconds!
    P.srate  = 32000;  % in sample pro second!
    P.figureHandle = [];
    P.figure = 1;
    P.binSize = 2;     % in milliseconds!
    P.axesHandles = [];
    P.figureTitle = [];
    P.subtractShiftPredictor = true;
    P.averageTrials = true;
    P.xc = [];
    P.xcs = [];
    P = mysort.util.parseInputs(P, 'xcorr', varargin);
    
    
    [spikeTrains IDs] = mysort.util.checkSpikeTrainsOrGDF(spiketrains_or_gdf, P);
%     spikeTrains = spikeTrains(~cellfun(@isempty, spikeTrains));
    if isempty(P.IDs) 
        P.IDs = IDs;
    end
    nSt     = length(P.IDs);
    nTrials = size(spikeTrains,1);
    %nXCorrs = nSt*(nSt-1)/2;
    assert(nSt>0, 'There must be at least one spikeTrain!');
    assert(nTrials>0, 'Some data would be nice, n''est pas?');
    assert(~P.subtractShiftPredictor || nTrials > 2, 'If the shift predictor has to be computed, more than 2 trials are necessary!');
    
    if isempty(P.axesHandles)
        if ~isempty(P.figureHandle)
            fh = P.figureHandle;
        elseif P.figure
            fh = mysort.plot.figure('name', 'XCorr');
        end  
        ah = mysort.plot.subplots(nSt*nSt,  'offsetY', .1, 'spacerY', 0, 'spacerX', 0); %'upperTriangle', 1,
    else
        ah = P.axesHandles;
        fh = get(ah(1), 'Parent');
    end
    
    if ~isempty(P.figureTitle)
        mysort.plot.figureTitle(P.figureTitle);
    end
    
    binSize_in_samples = P.binSize * P.srate/1000;
    maxLag_in_samples = P.maxLag * P.srate/1000;
    count = 1;
    for i=1:nSt
        id1 = P.IDs(i);
        for j=1:nSt            
            id2 = P.IDs(j);
            if P.subtractShiftPredictor
                if isempty(P.xc)
                    [xc xrange not_null_trials] = mysort.spiketrain.xcorr_raw(spikeTrains(1:end-1,i), spikeTrains(1:end-1,j),...
                        'maxLag', maxLag_in_samples, 'binSize', binSize_in_samples);
                else
                    xc = P.xc;
                end
                if isempty(P.xcs)
                    xcs = mysort.spiketrain.xcorr_raw(spikeTrains(1:end-1,i), spikeTrains(2:end,j),...
                        'maxLag', maxLag_in_samples, 'binSize', binSize_in_samples, 'T', P.T);
                else
                    xcs = P.xcs;
                end
                xc = xc-xcs;
            else
                if isempty(P.xc)
                    [xc xrange not_null_trials] = mysort.spiketrain.xcorr_raw(spikeTrains(:,i), spikeTrains(:,j),...
                        'maxLag', maxLag_in_samples, 'binSize', binSize_in_samples);
                else
                    xc = P.xc;
                end
            end
            nTrials = size(xc,1);
            if P.averageTrials
                xc = sum(xc,1)/nTrials;
            end
            
            xrange = xrange/(P.srate/1000);
            if ~isempty(xc)
                set(ah(count), 'NextPlot', 'add'); % that is hold on!
                plot(ah(count), xrange, xc, 'color', 'k');
            end
            if count == length(ah)
                xlabel(ah(count),'lag in ms');
            end
            if count == 1
                ylabel(ah(count),'xcorr');
            end
            axis tight

            if j==1
            else
                set(ah(count), 'xticklabel', []);
                set(ah(count), 'yticklabel', []);
            end
            
            if i==1
                title(ah(count), ['ID: ' num2str(id2)]);
                % Create ellipse
                pos = get(ah(count), 'position');
                pos(1) = pos(1) +  .85*pos(3);
                pos(2) = pos(2) + 1.01*pos(4);
                pos(3:4) = [.001 .001];
                elli = annotation(fh,'ellipse',pos,...
                         'FaceColor', mysort.plot.vectorColor(id2));
                set(elli, 'Units', 'pixel');
                pos = get(elli, 'position');
                set(elli, 'position', [pos(1:2) 10 10]);
                set(elli, 'Units', 'normalized');    
            end
            
            if j==nSt
                pos = get(ah(count), 'position');
                pos(1) = pos(1) + 1.01*pos(3);
                pos(2) = pos(2) +  .5 *pos(4);
                pos(3:4) = [.001 .001];
                elli = annotation(fh,'ellipse',pos,...
                         'FaceColor', mysort.plot.vectorColor(id1));
                set(elli, 'Units', 'pixel');
                pos = get(elli, 'position');
                set(elli, 'position', [pos(1:2) 10 10]);
                set(elli, 'Units', 'normalized');       
            end
            
            drawnow
            count = count +1;
        end
    end