function isiViolations(X, gdf, maxDist, varargin)
    P.maxSpikes = 1000;
    P.restrictTo = [];
    P.artefactEpochs = [];
    P.residualViolations = [];
    P.ax = [];
    P.fh = [];
    P = mysort.util.parseInputs(P, varargin);
    
    if ~isempty(P.restrictTo)
        gdf = gdf(ismember(gdf(:,1), P.restrictTo),:);
    end
    if isempty(gdf)
        return
    end
    units = unique(gdf(:,1));
    
    if isempty(P.ax)
        if isempty(P.fh)
            P.fh = mysort.plot.figure;
        end
        P.ax = axes();
    else
        P.fh = get(P.ax, 'parent');
    end
    set(P.ax, 'nextplot', 'add');
    for i=1:length(units)
        myts = gdf(gdf(:,1) == units(i),2);
        dts = diff(myts);
        violations = dts<maxDist;
        violations = [false; violations] | [violations; false];
        myts = myts(violations);
        plot(P.ax, myts, ones(length(myts),1)*i, 'o', 'markerfacecolor', mysort.plot.vectorColor(units(i)), 'markersize', 12);
    end
    if ~isempty(P.artefactEpochs)
        mysort.plot.epochs(P.artefactEpochs, -1, 'color', 'r', 'ax', P.ax);
    end
    if ~isempty(P.residualViolations)
        mysort.plot.epochs(P.residualViolations, -1.5, 'color', [0 0 .5], 'ax', P.ax);
    end
    set(P.ax, 'ylim', [-2 length(units)+1]);
    mysort.plot.dataCursorMode(P.fh);