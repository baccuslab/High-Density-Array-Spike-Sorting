function ax = axes(varargin)
    P.fontSize = 14;
    [P uP] = mysort.util.parseInputs(P, varargin, 'split');
    
    dP  = mysort.util.deflateP(P);
    duP = mysort.util.deflateP(uP);
    ax = axes(dP{:}, duP{:});