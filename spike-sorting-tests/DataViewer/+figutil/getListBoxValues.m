function [val, vals, idx, id] = getListBoxValues(handle, idlist)
    val = []; vals = []; idx = []; id = [];
    vals = get(handle,'String');
    idx  = get(handle,'Value');
    if iscell(vals) & ~isempty(idx) & idx(1) & ~isempty(vals) > 0
        val = vals(idx);
        if strcmp(val, '- - -')
            val = [];
        else
            if exist('idlist', 'var') && ~isempty(idlist)
                idx(idx>length(idlist)) = [];
                id = idlist(idx);
            end
        end     
    end