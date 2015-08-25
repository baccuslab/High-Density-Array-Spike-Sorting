function X = read_dset(dset_id, block, offset)
    %if size(offset,1)>size(offset,2); offset = offset'; end
    %if size(offset,1)>size(offset,2); offset = offset'; end
    assert(size(block,1) == 1, 'block must be a row vector!');
    assert(size(offset,1) == 1, 'offset must be a row vector!');
    
    nDims = length(block);
    assert(nDims == length(offset), 'Offset and number of block dims must match!');
    
    h5_offset = fliplr(offset); 
    h5_block  = fliplr(block); 
    h5_stride = [];
    h5_count  = ones(1, nDims);
    
    mem_space_id = H5S.create_simple(nDims, h5_block, h5_block);
    file_space_id = H5D.get_space(dset_id);      
    
    mode = 'H5S_SELECT_SET';
    plist = 'H5P_DEFAULT';
    H5S.select_hyperslab(file_space_id, mode, h5_offset, h5_stride, h5_count, h5_block);
    
    X = H5D.read(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, plist);