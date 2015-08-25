
function csvwrite(fname, M, append)
    str = [];
    if nargin < 3
        if exist(fname,'file'); delete(fname); end
    end
    for row=1:size(M,1)
        for col=1:size(M,2)
            if M(row,col) == round(M(row,col))
                str = [str sprintf('%d,', M(row,col))];
            else
                str = [str sprintf('%f,', M(row,col))];                
            end
        end
        str = [str(1:end-1) '\n'];
        if length(str) > 4000
           mysort.util.appendToFile(fname, str);
           str = [];
        end
    end
    
    mysort.util.appendToFile(fname, str);