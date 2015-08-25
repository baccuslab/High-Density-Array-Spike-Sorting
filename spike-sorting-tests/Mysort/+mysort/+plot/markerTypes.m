function m = markerTypes(i)
    m = {'o', 'x', 'd', 's'};
    if nargin == 0
        return
    end
    
    m = m{mod(i,length(m))+1};
        