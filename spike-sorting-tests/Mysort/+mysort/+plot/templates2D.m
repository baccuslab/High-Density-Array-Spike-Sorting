function P = templates2D(T, electrodePositions, maxChan, absThresh, varargin)
    P.fh = [];
    P.ah = [];
    P.stacked = 1;
    P.plotLegend = 1;
    P.sortTemplates = 'maxEl';
    P.IDs = [];
    P.nP = [];
    P = mysort.util.parseInputs(P, varargin);
    
    if isempty(P.ah) && isempty(P.fh)
        if P.stacked
            P.fh = mysort.plot.figure('w', 600, 'h', 900);
        else
            P.fh = mysort.plot.figure('w', 1100, 'h', 800);
        end
    end
    
    if isempty(P.nP)
        P.nP = size(T,3);
    end
    
    if isempty(P.ah) 
        if P.stacked
            P.ah = axes();
        else
            P.ah = mysort.plot.subplots(P.nP, 'preferVertical', 0);
        end
    end
    
    if isempty(P.IDs)
        P.IDs = 1:size(T,3);
    end
    
    if nargin < 3 || isempty(maxChan);
        maxChan = size(T,2);
    end
    if nargin < 4 || isempty(absThresh);
        absThresh = 0;
    end
    
    if strcmp(P.sortTemplates, 'maxEl')
        maxElPerTemplate = zeros(1, size(T,3));
        maximaPerTemplate = zeros( size(T,3), 1);
        for i=1:size(T,3)
            [maximaPerTemplate(i), maxElPerTemplate(i)] = max(max(abs(squeeze(T(:,:,i))), [], 1));
        end
        maxElPositionsPerTemplate = [electrodePositions(maxElPerTemplate,:) maximaPerTemplate];
        [B,idxmap] = sortrows(maxElPositionsPerTemplate, [2 1 -3]);
    else
        idxmap = 1:size(T,3);
    end
    legendstr = {};
    for i_=1:size(T,3)
        i = idxmap(i_);
        wfs = T(:,:,i);
        cidx = P.IDs(i);
        [col C marker] = mysort.plot.vectorColor(cidx);
        legendstr{i_} = num2str(cidx);
        if P.stacked==1
            P.LH{i} = mysort.plot.waveforms2D(wfs, electrodePositions, 'AxesHandle', P.ah, 'plotArgs', ...
                {marker, 'color', col, 'linewidth', 2},...
                'absThreshold', absThresh, 'maxNumberOfChannels', maxChan);
        else
            P.LH{i} = mysort.plot.waveforms2D(wfs, electrodePositions, 'AxesHandle', P.ah(i_), 'plotArgs', ...
                {'color', col, 'linewidth', 2},...
                'absThreshold', absThresh, 'maxNumberOfChannels', maxChan);            
            title(P.ah(i_), num2str(P.IDs(i)));
        end
    end
    if P.plotLegend
        legend(legendstr, 'location', 'EastOutside');
    end
    xlim = get(P.ah(1), 'xlim');
    ylim = get(P.ah(1), 'ylim');
    for i=2:min(length(P.ah), size(T,3))
        x = get(P.ah(i), 'xlim');
        y = get(P.ah(i), 'ylim');
        xlim = [min(xlim(1), x(1)) max(xlim(2), x(2))];
        ylim = [min(ylim(1), y(1)) max(ylim(2), y(2))];
    end
    set(P.ah, 'xlim', xlim, 'ylim', ylim);
    mysort.plot.figureName(P.fh, 'Templates 2D');
end