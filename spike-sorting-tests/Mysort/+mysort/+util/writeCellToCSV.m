
function writeCellToCSV(fname, C, append)
    str = [];
    if nargin < 3
        if exist(fname,'file'); delete(fname); end
    end
    for row=1:size(C,1)
        for col=1:size(C,2)
            if isnumeric(C{row,col})
                str = [str sprintf('%d,', C{row,col})];
            else
                str = [str sprintf('%s,', C{row,col})];                
            end
        end
        str = [str(1:end-1) '\n'];
    end
    
    mysort.util.appendToFile(fname, str);