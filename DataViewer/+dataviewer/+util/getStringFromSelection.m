function sel_str = getStringFromSelection(items)
    sel_str = '';
    if isempty(items)
        return
    end
    
    if iscell(items)
        for i=1:length(items)
            if ~isnumeric(items{i})
                items{i} = str2double(items{i});
            end
        end
        items = cell2mat(items);
    elseif ~isnumeric(items)
        sel_str = items;
        return
    end
    
    sel_str = util.dlmstring(items, ',');
    if length(sel_str) > 10
        a = num2str(min(items));
        b = num2str(max(items));
        sel_str = [a '...' b];
    end