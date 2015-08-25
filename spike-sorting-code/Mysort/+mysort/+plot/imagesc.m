% This is just a version of imagesc that makes nan values in X grey.
function [out cb] = imagesc(ah, X, varargin)
    if ~(max(size(ah)) == 1) || ~ishandle(ah)
        if nargin > 1
            varargin = [X varargin(:)'];
        end
        X = ah;
        ah = gca;
    end
    
    if ~any(isnan(X(:)))
        out = imagesc(X, 'parent', ah, varargin{:});
        cb_ = colorbar;
        return
    end
    
    map = colormap;
    map(1,:) = [.5 .5 .5];
    colormap(map);
    nColors = size(map,1);

    clims = [min(X(:)) max(X(:))];
    if clims(1)==clims(2);
        clims(1)=0;
    end
    X(isnan(X)) = clims(1) - (clims(2)-clims(1))/(nColors-2);
    out_ = imagesc(X, varargin{:});
    cb_ = colorbar('peer', ah);
    set(cb_, 'YLim' , clims);
    
    if nargout >= 1
        out = out_;
    end
    if nargout > 1
        cb = cb_;
    end
        