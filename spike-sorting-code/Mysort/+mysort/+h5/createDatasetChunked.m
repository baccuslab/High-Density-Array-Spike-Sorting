function datasetID = createDatasetChunked(fileID, dataspaceID, dsetname, h5type, chunkDims, deflation)
    datatypeID = H5T.copy(h5type); % H5T_NATIVE_FLOAT, H5T_NATIVE_DOUBLE, H5T_NATIVE_SHORT
    dcpl = H5P.create('H5P_DATASET_CREATE');
    lcpl = H5P.create('H5P_LINK_CREATE');
    H5P.set_create_intermediate_group(lcpl,1);
    dapl = 'H5P_DEFAULT';
    
    h5_chunk_dims = fliplr(chunkDims);
    H5P.set_chunk(dcpl, h5_chunk_dims);
    H5P.set_deflate(dcpl, deflation);
    datasetID = H5D.create(fileID, dsetname, datatypeID, dataspaceID, lcpl, dcpl, dapl);
    
    
%     datatypeID = H5T.copy(h5type); % H5T_NATIVE_FLOAT, H5T_NATIVE_DOUBLE, H5T_NATIVE_SHORT
%     dcpl = H5P.create('H5P_DATASET_CREATE'); % create property lis
%     h5_chunk_dims = fliplr(chunkDims);
%     H5P.set_chunk(dcpl, h5_chunk_dims);
%     H5P.set_deflate(dcpl, deflation);
%     datasetID = H5D.create(fileID, dsetname, datatypeID, dataspaceID, dcpl);
