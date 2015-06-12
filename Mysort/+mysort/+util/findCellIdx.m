function idx = findCellIdx(C, str)
    % finds the first index into the cell array C that contains the string
    % str
    idx = find(cellfun(@(x) strcmp(x, str), C),1);
end