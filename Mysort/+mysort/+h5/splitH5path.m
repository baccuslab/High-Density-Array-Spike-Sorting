function tok = splitH5path(h5path)
    expr = '/';
    [tok mat] = regexp(h5path, expr, 'split', 'match');
    for i=length(tok):-1:1
        if isempty(tok{i})
            tok(i) = [];
        end
    end