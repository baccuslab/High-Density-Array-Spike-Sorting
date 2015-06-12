
function waveforms(wavs, varargin)
	P.gdf = [];
	P.IDs = [];
    P.nC = 1;
    P.stacked = 1;
    P.axesHandle = [];
    P.plotMean = false;
    P.plotMedian = false;
    P.linewidth = 1;
    P.maxDist = [];
    P.restrict2Class = [];
    P.plotMaxNWaveformsPerClass = 1000;
	P = mysort.util.parseInputs(P, varargin, 'error');
    if isempty(wavs)
        return
    end
    if ndims(wavs) == 3
        [Tf P.nC nW] = size(wavs);
        wavs = mysort.wf.t2v(wavs);
    else
        nW = size(wavs,1);
        Tf = size(wavs,2)/P.nC;
    end

	
% 	assert(isempty(P.gdf) || isempty(P.IDs), 'Only gdf OR IDs may be provided!');
    
    if ~isempty(P.gdf)
        P.IDs = P.gdf(:,1);
    end
    if isempty(P.IDs)
        P.IDs = ones(nW,1);
    elseif length(P.IDs) == 1
        P.IDs = P.IDs*ones(nW,1);
    end
        
    if ~isempty(P.restrict2Class)
        idx = ismember(P.IDs, P.restrict2Class);
        wavs = wavs(idx,:);
        P.IDs = P.IDs(idx);
        if ~isempty(P.gdf)
            P.gdf = P.gdf(idx,:);
        end
    end
    
    if ~isempty(P.gdf)
        if ~isempty(P.maxDist);
            dt = [inf; diff(P.gdf(:,2))];
            idx = dt>P.maxDist;
            P.gdf(idx,:) = [];
            wavs(idx,:) = [];
            P.IDs(idx) = [];
        end
    end


    if ~P.stacked
        classes = unique(P.IDs);
        if isempty(P.axesHandle)
            fig = mysort.plot.figure();
            mysort.plot.figureName('Wfs');
            P.axesHandle = mysort.plot.subplots([length(classes) 1], 'spacerY', 0);
        end        
        glob_min = inf;
        glob_max = -inf;        
        for i=1:length(classes)
            glob_min = min(glob_min, min(min(wavs(P.IDs==classes(i),:))));
            glob_max = max(glob_max, max(max(wavs(P.IDs==classes(i),:))));
            mysort.plot.waveforms(wavs(P.IDs==classes(i),:), ...
                 P, 'gdf', [], 'linewidth', P.linewidth, ...
                'IDs', classes(i), 'stacked', 1,...
                'axesHandle', P.axesHandle(i), 'plotMaxNWaveformsPerClass', P.plotMaxNWaveformsPerClass);
            if i<length(classes)
                set(P.axesHandle(i), 'xtick', [], 'xticklabel', []);
            end
        end        
        set(P.axesHandle, 'ylim', [glob_min glob_max]);
        linkaxes(P.axesHandle, 'xy');
        return
    end
    
    if isempty(P.axesHandle)
        fig = mysort.plot.figure();
        mysort.plot.figureName('Wfs');
        P.axesHandle = axes();
    end
    
    wavs = mysort.wf.v4plot(wavs, P.nC);

    classes = unique(P.IDs);
    for i=1:length(classes)
        idx = find(P.IDs==classes(i));
        if length(idx) > P.plotMaxNWaveformsPerClass
            rp = randperm(length(idx));
            idx = idx(rp(1:P.plotMaxNWaveformsPerClass));
        end
        plot(P.axesHandle, wavs(idx,:)', ...
            'color', mysort.plot.vectorColor(classes(i)),...
            'linewidth', P.linewidth);
        set(P.axesHandle, 'nextplot', 'add');
        if P.plotMean
            plot(P.axesHandle, mean(wavs(P.IDs==classes(i),:), 1), 'color', 'k', 'linewidth', 3);
        end    
        if P.plotMedian
            plot(P.axesHandle, median(wavs(P.IDs==classes(i),:), 1), 'color', 'k', 'linewidth', 3);
        end                
    end
    
    mima = [min(wavs(:)) max(wavs(:))];
    for i=1:P.nC-1
        plot(P.axesHandle, [i*(Tf+1) i*(Tf+1)], mima, ':', 'color', [.6 .6 .8]);
    end
    axis(P.axesHandle, 'tight');
	