function b = isAlignedImage(X)
    b = isstruct(X) && ...
        isfield(X, 'xstart') && isfield(X, 'ystart') &&...
        isfield(X, 'xend') && isfield(X, 'yend') && isfield(X, 'IM');