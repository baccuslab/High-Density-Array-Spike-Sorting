
function templateSpikeTrain(T, gdf, varargin)
    P.mode = 'normal';
    P.cutLeft = 0;
    P.channelSpacer = 0;
    P.timeMultiplicator = 1;
    P.AxesHandle = [];
    P.srate = [];
    P.T_gdf_idx2id = [];  % use this variable if the ids in gdf(:,1) are not the indices into T
    [P unexplained_varargin] = mysort.util.parseInputs(P, varargin, 'split');
% 	[Tf nC nT] = size(T);
    
    if isempty(gdf)
        return
    end
    gdf = round(gdf);
    
    if isempty(P.AxesHandle)
        P.AxesHandle = axes();
        hold on
    end

    if ~isempty(P.srate)
        P.timeMultiplicator = 1/P.srate;
    end
    
    if strcmp(P.mode, 'normal')
        lw = 2;
        ls = '-';
        yoffset = 0;
    else %if strcmp(P.mode, 'groundtruth');
        lw = 3;
        ls = ':';
        yoffset = 0; %.1*abs(max(T(:)));
    end

    unexp = mysort.util.deflateP(unexplained_varargin);
    
%     id2idx = 1:size(T,3);
    idx2id = 1:size(T,3);
    if ~isempty(P.T_gdf_idx2id)
        idx2id = P.T_gdf_idx2id;
%         id2idx = @id2idx_helper;
    else
        ugdfid = unique(gdf(:,1));
        assert(~any(ugdfid > size(T,3)), 'If idx2id is not provided, the gdf ids must be indices into T !')
        P.T_gdf_id_idx = 1:size(T,3);
    end
    
    for idx=1:size(T,3)
        id = idx2id(idx);
        mygdf = gdf(gdf(:,1) == id,:);
        mygdf(:,1) = idx;
        if isempty(mygdf)
            continue
        end
        [x y] = mysort.wf.templateSpikeTrainData(T, mygdf, P.cutLeft); 
        if P.channelSpacer ~= 0
            for ch = 1:size(y,1)
                y(ch,:) = y(ch,:) + (ch-1)*P.channelSpacer;
            end
        end
        plot(P.AxesHandle, x*P.timeMultiplicator, yoffset+y', 'linewidth', lw, ...
                    'lineStyle', ls, ...
                    'color' , mysort.plot.vectorColor(id),...
                    unexp{:});
    end
%     %----------------------------------------------------------------------
%     function idx = id2idx_helper(id) 
%         idx = find(idx2id==id,1);
%     end
end