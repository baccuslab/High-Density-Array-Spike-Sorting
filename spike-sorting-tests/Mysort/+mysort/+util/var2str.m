
function str = var2str(var, varargin)
    P.formatForDouble = '%f';
    P.formatForInt = '%d';
    P = mysort.util.parseInputs(P, 'var2str', varargin);
    
    if iscell(var)
        str = var2str(var{1});
    elseif ischar(var)
        str = var;
    elseif var == int32(var)
        str = sprintf(P.formatForInt, var);
    else
        str = sprintf(P.formatForDouble, var);
    end
