function [sTemp bTemp selSpikes unselSpikes uTemp] = templatesGetSpikesForPlotting(T, C, t, nCplot)
    MAXPLOTSPIKES = 800;
    startTimeIdx = C.CutSpikesCutleft+1 - C.DisplaySpikesCutleft;
    %assert(startTimeIdx 
    selTimeIdx = startTimeIdx:startTimeIdx + C.DisplaySpikesTf -1;
    assert(selTimeIdx(1) >= 1, 'Cant display with the selected Cutleft/Tf combination (negative index)!');
    assert(selTimeIdx(end) <= C.CutSpikesTf, 'Cant display with the selected Cutleft/Tf combination (index is too big)!');
    
    nC = T.nC;
    if nargin <4
        nCplot = nC;
    end
    nCplot = min(nC, nCplot);
    
%     bTempM = mysort.wf.v2m(T.brushed_template(t,:), nC);
    sTempM = mysort.wf.v2m(T.getTemplateWf(t), nC);

    maxChans = max(abs(sTempM), [], 2);
    maxChans = sortrows([maxChans (1:nC)']);

    plotChanIdx = zeros(1,nC);
    plotChanIdx(maxChans(end-nCplot+1:end,2)) = 1;
    
    bTemp  = prune(T.getSourceTemplateWf(t), nC, plotChanIdx, selTimeIdx);
    sTemp  = prune(T.getTemplateWf(t), nC, plotChanIdx, selTimeIdx);
    uTemp  = prune(T.getExcludedTemplateWf(t), nC, plotChanIdx, selTimeIdx);    
    
    if nargout > 4
        t = T.getTemplateWithId(t);
        sIdx = t.getIdx;
        eIdx = t.getExcludedIdx();
        pIdx = randperm(t.getNWfs());
        pIdx = sIdx(pIdx(1:(min(MAXPLOTSPIKES, end))));
        selSpikes = prune(T.wfs.getWfs(pIdx), nC, plotChanIdx, selTimeIdx); 
        unselSpikes = [];
        if isempty(eIdx)
            return
        end
        pIdx = randperm(t.getNExcludedWfs());
        pIdx = eIdx(pIdx(1:(min(MAXPLOTSPIKES, end))));
        unselSpikes = prune(T.wfs.getWfs(pIdx), nC, plotChanIdx, selTimeIdx); 
    end
    %----------------------------------------------------------------------
    function V = prune(V, nC, cidx, tidx)
        nV = size(V,1);
        nCnew = sum(cidx);
        M = mysort.wf.v2m(V, nC);
        ccidx = repmat(cidx', nV, 1);
        V = mysort.wf.m2v([M(ccidx==1, tidx) nan(nV*nCnew,1)], nCnew);    
    end
end