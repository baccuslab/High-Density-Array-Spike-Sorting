function S = getAllChildrenWithTag(handle, type)
    % Retrieves all children of type "type" (optional) and stored the
    % handles as variables with the corresponding tag name in S
    S = [];
    if nargin == 2
        h = findobj(handle, 'type', type, '-regexp','Tag','[^'']');
        %findall(handle,'type','axes');
    else
        h = findobj(handle, '-regexp','Tag','[^'']');
        %h = get(handle, 'Children');
    end
    
    for i=1:length(h);
        S.(get(h(i), 'tag')) = h(i);
    end
   
