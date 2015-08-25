function logLastErrToFile(filename, path)

errStr = mysort.util.buildLastErrString();

if nargin < 2 path = []; end
[pathstr,name,ext] = fileparts(filename);

if strcmp(ext, '.h5') & ~isempty(path)
    % save to h5 file:
    l = length(errStr);    
    % number the errors:
    i = 1;
    while(i)
        numbered_path = [path num2str(i)];
        if ~mysort.h5.exist(filename, numbered_path)
            err_log = mysort.h5.createVariableAndOrFile(filename, numbered_path, [1 l], [1 l], 'H5T_C_S1');
            err_log(1,1:l) = errStr;
            clear err_log
            i = 0;
        end
    end
else
    % treat as simple text file:
    mysort.util.logToFile(filename, mysort.util.escapeString(errStr));
end

end