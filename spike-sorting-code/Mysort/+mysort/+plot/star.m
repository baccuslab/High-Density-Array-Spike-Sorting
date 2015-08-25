function star(X, cl, varargin)
    P.figureHandle = [];
    P.axesHandle = [];
    P.plotBasisVecs = 1;
    P.basisVecCols = [];
    P.markerForClass = {};
    P = mysort.util.parseInputs(P, 'plot.star', varargin);
    
    if isempty(P.axesHandle)
        if isempty(P.figureHandle)
            P.figureHandle = mysort.plot.figure('w', 800, 'h', 800);
            mysort.plot.figureName('StarPlot');
        end
        P.axesHandle = axes();
    end
    set(P.axesHandle, 'NextPlot', 'add');
    
    [S R order] = mysort.util.starTransform(X, 2);
    if ~exist('cl', 'var') || isempty(cl)
        cl = ones(size(X,1), 1);
    end
    
    if P.plotBasisVecs
        xlim = max(abs(S(:,1)));
        ylim = max(abs(S(:,2)));
        if isempty(P.basisVecCols)
            P.basisVecCols = repmat([.5 .5 .5], size(X,2), 1);
        end
        for i=1:size(X,2)
            if iscell(P.basisVecCols)
                col = P.basisVecCols{order(i)};
            else
                col = P.basisVecCols(order(i),:);
            end
            plot(P.axesHandle, [0 R(i,1)*xlim], [0 R(i,2).*ylim], 'color', [.5 .5 .5],...
                'linewidth', 2, 'color', col);
        end
    end
    
    cls = unique(cl);
    nCls = length(cls);
    for i=1:nCls
        idx = cl==cls(i);
        if ~isempty(P.markerForClass)
            mk = P.markerForClass{cls(i)};
        else
            mk = mysort.plot.markerForClass(cls(i));
        end
        plot(P.axesHandle, S(idx,1), S(idx,2), ... 
            'color', mysort.plot.vectorColor(i), ... 
            'marker', mk,...
            'linewidth', 2);
    end
    set(P.axesHandle, 'xlim', [-xlim xlim], 'ylim', [-ylim ylim]);

