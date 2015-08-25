function dataspaceID = createDataspace(dims, maxdims)
h5_dims = fliplr(dims);
nDims = length(dims);
if nargin == 1
    maxdims = dims;
else
    unlimited = H5ML.get_constant_value('H5S_UNLIMITED');
    maxdims(maxdims==-1) = unlimited;
end
h5_maxdims = fliplr(maxdims);

dataspaceID = H5S.create_simple(nDims, h5_dims, h5_maxdims); %{nC, 'H5S_UNLIMITED'}        