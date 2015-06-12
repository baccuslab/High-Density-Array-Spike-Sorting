
function v = toColVec(v)
    assert(size(v,1)==1 || size(v,2)==1, 'One dimension of v must have length 1!');
    assert(ndims(v) <=2, 'v must be a vector!');
    if size(v,1) < size(v,2)
        v = v';
    end
    