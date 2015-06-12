function [X yedges] = waveform22drep(y, yedges)
    nS = size(y,1);
    nX = size(y,2);
    
    if nargin < 2 || isempty(yedges)
        yedges = 50;
    end
    
    if length(yedges) == 1    
        yedges = linspace(min(y(:)), max(y(:)), yedges);
        yedges(end) = yedges(end)+1;
    end
    nY = length(yedges);
    
    X = zeros(nY, nX);
    
    % run over all x intervals
    for x = 1:nX-1
        X(:,x) = histc(y(:,x), yedges);    
    end