function X = h5read(fname, var_str, channelSet, idxSet)
    plist = 'H5P_DEFAULT';
    rmode = 'H5F_ACC_RDONLY';
    fid = H5F.open(fname, rmode, plist); 
    dset_id = H5D.open(fid, var_str);
    
    X = mysort.mea.h5read_dset(dset_id, channelSet, idxSet);
    
    H5D.close(dset_id);
    H5F.close(fid);