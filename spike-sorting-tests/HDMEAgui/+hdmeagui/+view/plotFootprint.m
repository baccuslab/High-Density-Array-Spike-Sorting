function plotFootprint(ax, T, t, C, GC)
    fprintf('  Plotting Footprint...\n');
    
    set(ax, 'NextPlot', 'replace');
    nC = T.nC;
    Tf = GC.DisplaySpikesTf;
    %% TODO !!! MULTI ELECTRODE
    if ~isa(T.wfs, 'mysort.wf.WfManagerBuffCut')
        return
    end
    CL = T.wfs.dataSource.getChannelList();
    CL = CL(CL(:,2)==1,3:4)/1000;

    
    [sTemp bTemp a b uTemp] = hdmeagui.view.templatesGetSpikesForPlotting(T, GC, t);
    bTempM = mysort.wf.v2m(bTemp, nC);
    sTempM = mysort.wf.v2m(sTemp, nC);
    uTempM = mysort.wf.v2m(uTemp, nC);
    
    % remove nans
    bTempM = bTempM(:, 1:end-1)/7;
    sTempM = sTempM(:, 1:end-1)/7;
    uTempM = uTempM(:, 1:end-1)/7;
    
    for i=1:nC
        idx = CL(i,1) + 15*(0:Tf-1)/Tf;
        if ~isempty(bTempM)
            y   = CL(i,2) - 1*uTempM(i,:);
            plot(ax, idx, y, '-', 'color', 'b', 'linewidth', 2);
            set(ax, 'NextPlot', 'add');
            y   = CL(i,2) - 1*bTempM(i,:);
            plot(ax, idx, y, '-', 'color', C.color_unselected_template, 'linewidth', 2);
        end
        set(ax, 'NextPlot', 'add');
        y   = CL(i,2) - 1*sTempM(i,:);
        plot(ax, idx, y, '-', 'color', C.color_selected_template, 'linewidth', 2);
    end

    axis(ax, 'tight');
    axis(ax, 'ij');
    fprintf('Update done.\n');
