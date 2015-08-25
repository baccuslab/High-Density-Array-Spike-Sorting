function mi_ma = minmax(x, mi_ma)
    if nargin == 1 || isempty(mi_ma)
        mi_ma = [min(x) max(x)];
    else
        mi_ma = [min(mi_ma(1),min(x)) max(mi_ma(2), max(x))];
    end