
function ah = xcorr(spiketrains_or_gdf, varargin)
    P.IDs = [];
    P.T = [];
    P.maxLag = 50;    % in milliseconds!
    P.srate  = 32000;  % in sample pro second!
    P.figureHandle = [];
    P.figure = 1;
    P.binSize = 2;     % in milliseconds!
    P.axesHandles = [];
    P.figureTitle = [];
    P.plotConfidence = true;
    P.subtractShiftPredictor = false; % Dont use this
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
    maxAxis = -inf;
    minAxis =  inf;
    for i=1:nSt
        id1 = P.IDs(i);
        for j=1:nSt            
            id2 = P.IDs(j);
            [xcorr E_xc var_xc conf_95_xc xrange edg xc_unnormalized] = mysort.spiketrain.xcorr(spikeTrains(:,i), spikeTrains(:,j),...
                'maxLag', maxLag_in_samples, 'binSize', binSize_in_samples, 'T', P.T);
            if P.subtractShiftPredictor
                error('nor properly implemented. Use xcorr_raw to do this');
                [xcorrs E_xcs var_xcs conf_95_xcs xranges edgs xc_unnormalizeds] = mysort.spiketrain.xcorr(spikeTrains(1:end-1,i), spikeTrains(2:end,j),...
                    'maxLag', maxLag_in_samples, 'binSize', binSize_in_samples, 'T', P.T);
                xcorr = xcorr-xcorrs;
            end
%             [XCorr]=trains.brillingerXCorr( spikeTrains{i}(:), spikeTrains{j}(:), binSize_in_samples, maxLag_in_samples);
            
%             figure;
%             subplot(3,2,1)
%             plot(xc_unnormalized); hold on
%             subplot(3,2,2)
%             plot(xcorr); hold on
%             subplot(3,2,3)
%             plot(XCorr.p_est_mn, 'g')
%             subplot(3,2,4)
%             plot(XCorr.normed, 'g')
            
            xrange = xrange/(P.srate/1000);
            if ~isempty(xcorr)
                xcorr  = sqrt(xcorr);
                E_xc   = sqrt(E_xc);
                var_xc = sqrt(var_xc);
                conf_95_xc = sqrt(conf_95_xc);

                xcorr  = xcorr - E_xc;
                old_E_xc = - E_xc;
                E_xc   = 0;
                xcorr  = xcorr ./ var_xc;
                conf_95_xc = conf_95_xc./ var_xc;
                old_var_xc = var_xc;
                var_xc = 1;
                
                plot(ah(count), xrange([1 end]), E_xc*[1 1], 'color', 'k', 'linestyle', ':');
                set(ah(count), 'NextPlot', 'add'); % that is hold on!
                plot(ah(count), xrange, xcorr, 'color', 'k');
                if P.plotConfidence
                    plot(ah(count), xrange([1 end]), E_xc*[1 1]+conf_95_xc, 'color', 'r', 'linestyle', ':');
                    plot(ah(count), xrange([1 end]), E_xc*[1 1]-conf_95_xc, 'color', 'r', 'linestyle', ':');
                end
                plot(ah(count), xrange([1 end]), old_E_xc/old_var_xc*[1 1], 'color', 'b', 'linestyle', '-');
            end
            if count == length(ah)
                xlabel(ah(count),'lag in ms');
            end
            if count == 1
                ylabel(ah(count),'xcorr');
            end
            axis tight

            if ~isempty(xcorr)
                if i==j
                    valIdx = [1:(length(xcorr)-1)/2 2+(length(xcorr)-1)/2:length(xcorr)];
                    mi_ma = [min(xcorr(valIdx))  max(xcorr(valIdx))];
                    if mi_ma(1) == mi_ma(2); mi_ma(2) = xcorr((length(xcorr)-1)/2 +1 ); end
                    if mi_ma(1) == mi_ma(2); mi_ma = [0 1]; end
                    set(ah(count), 'ylim', mi_ma);
                elseif ~isempty(xcorr)
                    mi_ma = [min(xcorr) max(xcorr)];
                end
                minAxis = min(minAxis, mi_ma(1));
                maxAxis = max(maxAxis, mi_ma(2));                
            end
            
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
    if ~isinf(minAxis)
        set(ah, 'ylim', [minAxis maxAxis]);
    end