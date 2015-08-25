
function writeCell(fname, tree, c, mode)
    if nargin == 3
        mode = 'overwrite';
    end
    
    for i=1:length(c)
        mysort.h5.recursiveSave(fname, c{i}, [tree '/' sprintf('%05d',i)], mode);
    end
    