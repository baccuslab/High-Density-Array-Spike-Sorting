function X = write_hyperslab(fname, h5path, X, offset)
    plist = 'H5P_DEFAULT';
    rmode = 'H5F_ACC_RDWR';
    fid = H5F.open(fname, rmode, plist); 
    dset_id = H5D.open(fid, h5path);
    mysort.h5.write_dset(dset_id, X, offset);
    H5D.close(dset_id);
    H5F.close(fid);
    
