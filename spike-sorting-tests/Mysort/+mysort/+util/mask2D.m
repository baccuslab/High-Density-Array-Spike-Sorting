function neighbors = mask2D(X, my_fun, neighbor_fun, neighborhood, minNeighbors)
    nX = size(X,2);
    nY = size(X,1);
    
    if nargin < 4 || isempty(neighborhood)
        neighborhood = [-1 -1
               -1  0
               -1  1
               0  -1
               0  1
               1  -1
               1  0
               1  1];
    end
    neighbors = zeros(size(X));
    for xidx = 1:nX
        for yidx = 1:nY
            if my_fun(xidx, yidx)
                for k=1:size(neighborhood,1)
                    Nx = neighborhood(k,1);
                    Ny = neighborhood(k,2);
                    if xidx + Nx > 0   && yidx + Ny > 0   && ...
                       xidx + Nx <= nX && yidx + Ny <= nY && ...   
                       neighbor_fun(xidx+Nx, yidx+Ny)
                        neighbors(xidx, yidx) = neighbors(xidx, yidx) + 1;
                    end
                end
            end
        end
    end
    % mask the mask again, so that all pixels that actually contributed to
    % one pixel being over threshold are also kept. this avoids shrinking
    % of the mask
    if nargin > 4 && ~isempty(minNeighbors)
        neighbors = mysort.util.mask2D(X, my_fun, neighbors>=minNeighbors, neighborhood);
    end