
function epochs(epochs, y, varargin)
    P.ax = [];
    P.fh = [];
    P.linewidth = 8;
    [P uP] = mysort.util.parseInputs(P, varargin, 'split');
    uP = mysort.util.deflateP(uP);
    
    if isempty(P.fh) && isempty(P.ax)
        fh = mysort.plot.figure;
        P.ax = axes();
    elseif isempty(P.ax)
        isempty(P.ax)
        P.ax = axes();
    end
        
    
    if ~exist('y', 'var')
        y = 1;
    end

    o = y*[ones(size(epochs)) nan(size(epochs,1),1)];
    epochs = [epochs nan(size(epochs,1), 1)];

    epochs = epochs';
    epochs = epochs(:);
    o = o';
    o = o(:);
    
    plot(P.ax, epochs, o, uP{:}, 'linewidth', P.linewidth);