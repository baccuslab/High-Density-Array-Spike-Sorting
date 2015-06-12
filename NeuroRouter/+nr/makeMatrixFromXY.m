function M = makeMatrixFromXY(x,y) %, varargin)
%     P.xtolerance = 3;
%     P.ytolerance = 3;
%     P = mysort.util.parseInputs(P, varargin, 'error');
    
    nP = length(x);
    assert(nP == length(y), 'x and y must have same lengths')
    
    % first, identify unique groups of x positions and y positions inside
    % tolerane
    
%     x = round(x/P.xtolerance);
    ux = unique(x);
    
%     y = round(y/P.ytolerance);
    uy = unique(y);
    
    nX = length(ux);
    nY = length(uy);
    
    M = zeros(nX,nY);
    
    for i=1:nP
        xidx = find(ux==x(i),1);
        yidx = find(uy==y(i),1);
        M(xidx,yidx) = i;
    end
    
    
    
    