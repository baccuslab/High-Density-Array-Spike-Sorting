function b = dataSourceIsEmpty(ds)
    b = true;
    if isempty(ds)
        return
    end
    if iscell(ds)
        b = isempty(ds{1});
    end
        