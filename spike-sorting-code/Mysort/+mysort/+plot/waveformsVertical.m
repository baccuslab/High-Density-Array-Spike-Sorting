
function P = waveformsVertical(wavs, varargin)
	P.gdf = [];
	P.IDs = [];
    P.nC = [];
    P.stacked = 1;
    P.axesHandle = [];
    P.plotMean = false;
    P.plotMedian = false;
    P.linewidth = [];
    P.channelSpacer = [];
    P.channelSpaceFactor = .75;
    P.maxDist = [];
    P.restrict2Class = [];
    P.plotMaxNWaveformsPerClass = 1000;
    P.xOffset = 0;
    P.forceV = 0;
    P.zeroLineMarker = 'k:';
    P.ignoreClassesForColoring = 0;
    P.plotArgs = {'-', 'color', [.5 .5 .5]};
    
	P = mysort.util.parseInputs(P, varargin, 'error');
    if isempty(wavs)
        return
    end
    % for backward compatibility
    if ~isempty(P.linewidth)
        P.plotArgs = [P.plotArgs {'linewidth', P.linewidth}];
    end
    
    wfsIsConcatenated = P.forceV==1 || (ismatrix(wavs) && ~isempty(P.nC));
    if ~wfsIsConcatenated
        if ndims(wavs) < 3
            warning('This function only accepts wfs in tensor form. Use forceV switch to use wfs in v form.');
        end
        [Tf P.nC nW] = size(wavs);
        wavs = mysort.wf.t2v(wavs);
    else
        assert(~isempty(P.nC), 'If wavs are concatenated waveforms the number of channels must be provided in P.nC!');
        nW = size(wavs,1);
        Tf = size(wavs,2)/P.nC;
    end

	
% 	assert(isempty(P.gdf) || isempty(P.IDs), 'Only gdf OR IDs may be provided!');

    if ~iscell(P.plotArgs)
        P.plotArgs = {P.plotArgs};
    end

    if ~isempty(P.gdf)
        P.IDs = P.gdf(:,1);
    end
    if isempty(P.IDs)
        P.IDs = ones(nW,1);
    elseif length(P.IDs) == 1
        P.IDs = P.IDs*ones(nW,1);
    end
    assert(length(P.IDs) == nW, 'Not the right number of IDs!');
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
    min_maxis = zeros(2, P.nC);
    channel_spacer = zeros(1, P.nC);
    for c=1:P.nC
        cidx = mysort.wf.vSubChannelIdx(wavs, P.nC, c);
        m1 = max(max(wavs(:,cidx)));
        m2 = min(min(wavs(:,cidx)));
        min_maxis(:,c) = [m1; m2];
        if isempty(P.channelSpacer)
            if c>1
                channel_spacer(c) = P.channelSpaceFactor*(min_maxis(1, c-1) - min_maxis(2, c));
            end
        elseif length(P.channelSpacer) == 1
            channel_spacer(c) = P.channelSpacer;
        else
            channel_spacer(c) = P.channelSpacer(c);
        end        
    end

    if ~P.stacked
        classes = unique(P.IDs);
        if isempty(P.axesHandle)
            fig = mysort.plot.figure();
            mysort.plot.figureName('Wfs');
            P.axesHandle = mysort.plot.subplots([1 length(classes)], 'spacerY', 0);
        end        
        for i=1:length(classes)
            mysort.plot.waveformsVertical(wavs(P.IDs==classes(i),:), ...
                 P, 'gdf', [], ...
                'IDs', classes(i), 'stacked', 1,...
                'channelSpacer', channel_spacer,...
                'xOffset', P.xOffset,...
                'forceV', 1,...
                'plotArgs', P.plotArgs,...
                'axesHandle', P.axesHandle(i), 'plotMaxNWaveformsPerClass', P.plotMaxNWaveformsPerClass);
%             if i<length(classes)
%                 set(P.axesHandle(i), 'xtick', [], 'xticklabel', []);
%             end
            if i>1
                set(P.axesHandle(i), 'ytick', [], 'yticklabel', []);
            end
        end        
        set(P.axesHandle, 'ylim', [min_maxis(2,1) sum(channel_spacer)+min_maxis(1,end)]);
        linkaxes(P.axesHandle, 'xy');
        return
    end
    
    if isempty(P.axesHandle)
        fig = mysort.plot.figure();
        mysort.plot.figureName('Wfs');
        P.axesHandle = axes();
    end
    for c=1:P.nC
        cidx = mysort.wf.vSubChannelIdx(wavs, P.nC, c);
        wavs(:,cidx) = wavs(:,cidx) + sum(channel_spacer(1:c));
    end    
    wavs = mysort.wf.v4plot(wavs, P.nC);
    classes = unique(P.IDs);
    xidx = repmat([1:Tf nan], 1, P.nC)+P.xOffset;

    for i=1:length(classes)
        idx = find(P.IDs==classes(i));
        if length(idx) > P.plotMaxNWaveformsPerClass
            rp = randperm(length(idx));
            idx = idx(rp(1:P.plotMaxNWaveformsPerClass));
        end
        
        if P.ignoreClassesForColoring
            plot(P.axesHandle, xidx, wavs(idx,:)', ...
                P.plotArgs{:});
        else
            plot(P.axesHandle, xidx, wavs(idx,:)', ...
                P.plotArgs{:}, 'color', mysort.plot.vectorColor(classes(i)));
        end
        set(P.axesHandle, 'nextplot', 'add');
        if ~isempty(P.zeroLineMarker)
            for c=1:P.nC
                plot(P.axesHandle, [1 Tf], ones(1,2)*sum(channel_spacer(1:c)), P.zeroLineMarker);
            end        
        end
%         if P.plotMean
%             plot(P.axesHandle, mean(wavs(P.IDs==classes(i),:), 1), 'color', 'k', 'linewidth', 3);
%         end    
        if P.plotMedian
            plot(P.axesHandle, xidx, median(wavs(P.IDs==classes(i),:), 1), 'color', 'k', 'linewidth', 3);
        end             
    end
    axis(P.axesHandle, 'tight');
    P.channelSpacer = channel_spacer;
	