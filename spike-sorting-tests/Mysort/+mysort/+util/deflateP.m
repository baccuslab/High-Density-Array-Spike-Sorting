
function args = deflateP(P)
    if isempty(P)
        args = {};
        return
    end
    if ~isstruct(P)
        args = P;
        return
    end
    fnames = fieldnames(P);
    args = cell(1, length(fnames));
    for i=1:length(fnames);
        args{2*i -1} = fnames{i};
        args{2*i   } = P.(fnames{i});
    end