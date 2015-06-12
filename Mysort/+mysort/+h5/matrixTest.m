clear M
clear mysort.h5.matrix

fname = 'matrixTest4.h5';
h5path = '/foo/bar';

if exist(fname, 'file')
    M = mysort.h5.matrix(fname, h5path);
else
    dims = [1000 10];
    maxDims = [H5ML.get_constant_value('H5S_UNLIMITED') 10];
    h5type = 'H5T_NATIVE_DOUBLE';
    M = mysort.h5.matrixCreate(fname, h5path, dims, maxDims, h5type);   
end

M(1,1)
M(1,1) = -1;
M(1,1)
