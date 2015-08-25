
function figureName(fh, str)
    if nargin == 1
        str = fh;
        fh = gcf;
    end
    set(gcf, 'Name', str);