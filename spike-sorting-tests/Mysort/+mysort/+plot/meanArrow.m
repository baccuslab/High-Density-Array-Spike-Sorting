function P = meanArrow(M, varargin)
    P.ah = [];
    P = mysort.util.parseInputs(P, varargin);
    
    
    if isempty(P.ah)
        P.ah = gca;
    end
    set(P.ah, 'nextplot', 'add');
    lims = get(P.ah, 'ylim');
    plot(P.ah, [M M], [.95*lims(2) lims(2)], 'k-', 'linewidth', 2);
    

    