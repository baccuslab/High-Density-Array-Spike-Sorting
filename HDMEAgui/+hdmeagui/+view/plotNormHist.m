function plotNormHist(ax, G, t, C, GC)    
    T = G.T;
    fprintf('  Plotting NormHist...\n');
    set(ax, 'NextPlot', 'replacechildren');
    P = G.getCurrentTemplateProjections();
    Pex = G.getCurrentTemplateExcludedProjections();
    if isempty(Pex) 
        e = mysort.plot.normalityHist(P,...
            'axesHandles', ax, ...
            'normFitPara',{'color', 'r', 'linewidth', 2},...
            'plotNormFit', 0, 'nBins', 40);
    else
        e = mysort.plot.normalityHist([Pex; P],...
            'axesHandles', ax, 'nBins', 40, ...
            'plotNormFit', 0, 'facecolor', C.color_unselected_spikes);
        set(ax, 'NextPlot', 'add');
        mysort.plot.normalityHist(P,...
            'axesHandles', ax, 'nBins', 40, ...
            'normFitPara',{'color', C.color_selected_template, 'linewidth', 2},...
            'plotNormFit', 0, 'facecolor', C.color_selected_spikes, 'edges', e);
    end
    set(ax, 'xlim', [e(1)*.9 e(end)])
    set(ax, 'NextPlot', 'add');
    axes(ax);
    G.normHistPlotThresholdHandle = ...
        mysort.plot.verticalLines(G.normHistThresh, [], C.color_normhist_old_thresh);
    ch = get(ax, 'children');
    set(ch, 'Hittest', 'off');
    xlabel(ax, 'projection on template');
    ylabel(ax, 'count');
    