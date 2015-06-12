function X = read_dset_old(dset_id, channelSet, idxSet)
    c1 = min(channelSet);
    c2 = max(channelSet);
    i1 = min(idxSet);
    i2 = max(idxSet);
    nC_ = c2-c1+1;
    nS_ = i2-i1+1;
    assert(nC_*nS_ < 10^9, sprintf('Loading that much data (%d) at once is not recommended!', nC_*nS_));
    block = [nC_ nS_];
    h5_block = block; 
    mode = 'H5S_SELECT_SET';
    plist = 'H5P_DEFAULT';
    % read outer bounding box of requested data, this is usually faster
    % than reading individual rows or indices
    dims =  [nC_ nS_];  
    h5_dims = dims; 
    mem_space_id = H5S.create_simple(2, h5_dims, h5_dims);
    file_space_id = H5D.get_space(dset_id);
    offset = [c1-1 i1-1];      
    h5_offset = offset; 
    H5S.select_hyperslab(file_space_id,mode,h5_offset,[],[1 1],h5_block);
    X = H5D.read(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, plist);

    % select actually requested data
    if length(idxSet) < nS_ && length(channelSet) < nC_
        X = X(idxSet-i1+1, channelSet-c1+1);
    elseif length(idxSet) < nS_ 
        X = X(idxSet-i1+1, :);
    elseif length(channelSet) < nC_
        X = X(:, channelSet-c1+1);
    end
