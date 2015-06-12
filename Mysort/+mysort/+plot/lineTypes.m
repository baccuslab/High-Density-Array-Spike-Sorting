function m = lineTypes(i)
    m = {'-', ':', '--', '-.'};
    if nargin == 0
        return
    end
    
    m = m{mod(i,length(m))+1};
        