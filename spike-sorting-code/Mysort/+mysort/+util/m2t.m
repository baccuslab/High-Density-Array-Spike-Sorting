
function T = m2t(M, nC)
    warning('This function is depricated! Use mysort.wf.* instead!');
    [items, dim] = size(M);
    dim = dim/nC;
    assert(round(dim)==dim, 'Number of channels does not match!');
    T = zeros(dim,nC,items);
    for i=1:items
        T(:,:,i) = reshape(M(i,:)', dim, nC);
    end
