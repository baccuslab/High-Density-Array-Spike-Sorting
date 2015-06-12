function P = vectors(V, varargin)
    P.fh = [];
    P.ah = [];
    P.plotArgs = {};
    P.autoColor = 1;
    P = mysort.util.parseInputs(P, varargin, 'error');
    if isempty(P.ah)
        if isempty(P.fh)
            P.fh = mysort.plot.figure([500 500]);
        end
        P.ah = axes(); 
        
    end
    plot(P.ah, 0,0, 'k+', 'markersize', 18, 'linewidth', 2)
    hold on
    if P.autoColor
        for i=1:size(V,1)
            plot(P.ah, [0 V(i,1)], [0 V(i,2)], P.plotArgs{:}, 'color', mysort.plot.vectorColor(i));
        end
    else
        for i=1:size(V,1)
            plot(P.ah, [0 V(i,1)], [0 V(i,2)], P.plotArgs{:});
        end
    end
    axis(P.ah, 'equal')