
function GDFwrite(fileName, X, varargin)
    % Writes an Array to a .gdf File
    % Inputs: 
    %   fileName - Path and Name to .gdf File, will be created if
    %              nonexistent overwritten otherwise
    %   X        - Datamatrix, first Column ID, second Timestamp
    P.delimiter = 'space';
    P = mysort.util.parseInputs(P,'GDFwrite',varargin);
    
    X = X';
    delimiter = ' ';
    if strcmp(P.delimiter, 'tab')
        delimiter = sprintf('\t');
    end
    fhandle = fopen(fileName, 'w+');
    fprintf(fhandle, ['%05d' delimiter '%d\n'], X(:));
    fclose(fhandle);