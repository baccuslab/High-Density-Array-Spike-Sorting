function datasetID = createDataset(fileID, dataspaceID, dsetname, h5type)
    datatypeID = H5T.copy(h5type); % H5T_NATIVE_FLOAT, H5T_NATIVE_DOUBLE, H5T_NATIVE_SHORT, H5T_C_S1
%     if strcmp(h5type, 'H5T_C_S1')
%         % This might be necessary to have the H5 viewer depict the string
%         % correctly
%         H5T.set_size(datatypeID,'H5T_VARIABLE');
%     end
%     
    dcpl = H5P.create('H5P_DATASET_CREATE');
    lcpl = H5P.create('H5P_LINK_CREATE');
    H5P.set_create_intermediate_group(lcpl,1);
    dapl = 'H5P_DEFAULT';

    datasetID = H5D.create(fileID, dsetname, datatypeID, dataspaceID, lcpl, dcpl, dapl);
