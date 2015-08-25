function [b subH5info Dims bIsArray] = exist(filename, h5path, h5infoStruct)
    % This function checks if inside an existing h5 file a dataset with the
    % specified path exists
    bIsArray = false;
    Dims = [];
    try 
        if ~isempty(strfind(version, 'R2012'))
            subH5info = h5info(filename, h5path);
            b=true;
            bIsArray = isfield(subH5info, 'Dataspace') && max(subH5info.Dataspace.Size) > 1;
            Dims = subH5info.Dataspace.Size;
            return
        else
            if ~exist(filename, 'file')
                b = false;
                subH5info = [];
                return
            end
            if nargin < 3 || isempty(h5infoStruct)
                h5infoStruct = hdf5info(filename);
            end
            subH5info = mysort.h5.findSubH5info(h5infoStruct, h5path);
            b = ~isempty(subH5info);            
            if ~isempty(subH5info)
                bIsArray = isfield(subH5info, 'Dims') && max(subH5info.Dims) > 1;
                Dims = subH5info.Dims;
            end
        end
    catch
        subH5info = [];
        b=false;
    end
