function idx = setUserSelection(handle, string_values)
    if ~iscell(string_values)
        string_values = {string_values};
    end
    possible_vals = get(handle, 'String');
    idx = [];
    for i=1:size(possible_vals,1)
        for j=1:length(string_values)
            if strcmp(possible_vals{i}, string_values{j})
                idx = [idx i];
            end
        end
    end        
    if ~isempty(idx)
        set(handle, 'Value', idx);
    end

          