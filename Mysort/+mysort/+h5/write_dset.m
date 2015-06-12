function write_dset(dset_id, X, offset)
    if size(offset,1)>size(offset,2); offset = offset'; end
    blockDims = size(X);
    if length(offset) == 1
        assert(any(blockDims==1), 'One Blockdim must have size 1 if 1 dim array!');
        blockDims = max(blockDims);
    end
    nDims = length(blockDims);
    assert(nDims==size(offset,2), 'Offset and number of block dims must match!');
    assert(size(offset,1) == 1, 'Offset must be a vector of length nDims!');
    h5_start = fliplr(offset);
    h5_block = fliplr(blockDims);
    h5_stride = [];
    h5_count  = ones(1, nDims);
    
    mem_space_id = H5S.create_simple(nDims, h5_block, []);
    file_space_id = H5D.get_space(dset_id);

    H5S.select_hyperslab(file_space_id, 'H5S_SELECT_SET', h5_start, ...
        h5_stride, h5_count, h5_block);

    H5D.write(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, 'H5P_DEFAULT', X);
    
