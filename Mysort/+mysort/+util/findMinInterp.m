function [x fval idx m] = findMinInterp(y)
    % finds the interpolated minimum of y. uses a heuristic to find the
    % interval in y that could contain the minimum and uses fminbnd on an
    % interpolation function handle to actually find it
    
    [m idx] = min(y);
    
    interp_range    = max(1,idx-5):min(length(y), idx+5);
    search_range_start = max(1,idx-1) - interp_range(1)+1;
    search_range_end   = min(length(y), idx+1) - interp_range(1) +1;
    
    opt = [];
    opt.MaxFunEvals = 30;
    opt.TolX = 0.01;
    opt.Display = 'off';    
    
    ipf = mysort.util.interpolfun(y(interp_range));
    [x, fval] = fminbnd(ipf, search_range_start, search_range_end, opt);
    x = x + interp_range(1)-1;
    
    
    
    
    