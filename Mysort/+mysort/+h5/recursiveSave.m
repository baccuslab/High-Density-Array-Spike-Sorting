function recursiveSave(fname, val, tree, mode)
    if nargin == 2
        tree = [];
        mode = 'overwrite';
    end
    if nargin == 3
        mode = 'overwrite';
    end    
    if strcmp(mode, 'overwrite')
        hdf5write(fname, '/', [], 'WriteMode', 'overwrite');
    end    
    varMode = 'append';
    if isstruct(val)
        flist = fieldnames(val);
        if length(val) > 1
            % Struct Array!
            for k=1:length(val)
                myname = sprintf('%05d/', k);
                for i=1:length(flist)
                    mysort.h5.recursiveSave(fname, val(k).(flist{i}), [tree '/' myname flist{i}], varMode);
                end
            end
        else
            for i=1:length(flist)
                mysort.h5.recursiveSave(fname, val.(flist{i}), [tree '/' flist{i}], varMode);
            end
        end
    else
        if iscell(val)
            mysort.h5.writeCell(fname, tree, val, varMode);
        elseif isa(val, 'mysort.util.DebuggableClass')
            mysort.h5.recursiveSave(fname, val.getStruct(), [tree '____CLASS____' class(val)], varMode);
        elseif islogical(val)
            check_create()
            hdf5write(fname, tree, double(val), 'WriteMode', varMode, 'V71Dimensions', true);
        elseif isa(val, 'function_handle')
            check_create()
            hdf5write(fname, tree, func2str(val), 'WriteMode', varMode, 'V71Dimensions', true);
        elseif isempty(val)
            check_create()
            hdf5write(fname, [tree '____EMPTY____'], '____EMPTY____', 'WriteMode', varMode, 'V71Dimensions', true); 
        else
            check_create()
            hdf5write(fname, tree, val, 'WriteMode', varMode, 'V71Dimensions', true);
        end
    end
    
    function check_create()
        if ~exist(fname, 'file')
            hdf5write(fname, '/createdByRecursiveSave', 1, 'WriteMode', 'overwrite');
        end
    end
end