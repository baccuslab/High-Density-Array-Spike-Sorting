function M = createVariableAndOrFile(fname, h5path, dims, maxDims, h5type, chunkDims, deflation, bReadOnly)
    if nargin < 8
        bReadOnly = false;
    end
    if nargin == 5 || isempty(chunkDims)
        unlimited = H5ML.get_constant_value('H5S_UNLIMITED'); 
        if any(maxDims == unlimited)
            warning('specify chunk dimensions for unlimited datasets!\nTaking default');
            chunkDims = maxDims;
            chunkDims(maxDims==unlimited) = 1000;
            deflation = 0;
        else
            chunkDims = [];
            deflation = [];
        end
    end
    % H5T_NATIVE_FLOAT, H5T_NATIVE_DOUBLE, H5T_NATIVE_SHORT

    if exist(fname, 'file')
%         disp('File already exist, trying to create h5 variable...')
        plist = 'H5P_DEFAULT';
        if bReadOnly
            rmode = 'H5F_ACC_RDONLY';
        else
            rmode = 'H5F_ACC_RDWR';
        end
        fileID = H5F.open(fname, rmode, plist); 
    else
        disp('Creating File...')
        fileID = mysort.h5.createFile(fname);
    end
    dataspaceID = mysort.h5.createDataspace(dims, maxDims);
    
    if isempty(chunkDims)
        datasetID = mysort.h5.createDataset(fileID, ...
            dataspaceID, h5path, h5type);
    else
        datasetID = mysort.h5.createDatasetChunked(fileID, ...
            dataspaceID, h5path, h5type, ...
            chunkDims, deflation);
    end
    H5D.close(datasetID);
    H5F.close(fileID); 
    M = mysort.h5.matrix(fname, h5path, bReadOnly);