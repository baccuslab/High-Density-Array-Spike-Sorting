
function h = plotArray(x, y, varargin)
    P.axesHandle = [];
    P.cla = 1;
    P.color = 'k';
    [P uP] = mysort.util.parseInputs(P, varargin, 'split');
    
    if isempty(P.axesHandle)
        P.axesHandle = axes();
    end
    if P.cla
        cla(P.axesHandle);
    end
    set(P.axesHandle, 'nextPlot', 'add');
    tmp = mysort.util.deflateP(uP);
    h = plot(P.axesHandle, x, y, '.', 'color', P.color, tmp{:});
    
    % axis ij
    set(P.axesHandle,...
        'XDir','normal',...
        'YDir','reverse');