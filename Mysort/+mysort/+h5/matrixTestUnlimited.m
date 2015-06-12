clear M
clear mysort.h5.matrix

fname = 'matrixTestUnlimited.h5';
h5path = 'foo'; %/foo/bar';

if exist(fname, 'file')
    M = mysort.h5.matrix(fname, h5path);
else
    unlimited = H5ML.get_constant_value('H5S_UNLIMITED');
    dims = [10 10];
    maxDims = [-1 -1];
    h5type = 'H5T_NATIVE_DOUBLE';
    M = mysort.h5.matrixCreate(fname, h5path, dims, maxDims, h5type);   
end

M(1,1)
M(1,1) = -1;
M(1,1)


M(11,11) = -1;
M(:,:)
