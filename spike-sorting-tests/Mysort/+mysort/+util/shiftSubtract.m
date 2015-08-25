
function x = shiftSubtract(x, y, offset, inf_flag)
    % Subtracts the possible multi channel vector (rows are channels) y 
    % from the vector x with a certain offset. x and y can due to the 
    % offset be only partly overlapping
    if isempty(x)
        return
    end
    if exist('inf_flag', 'var') && inf_flag==1
        func = @inf;
    elseif exist('inf_flag', 'var') && inf_flag==-1
        func = @(a,b) -inf(a,b);
    elseif exist('inf_flag', 'var') && inf_flag==2
        func = @nan;
    else
        func = @zeros;
    end
    assert(~isempty(y), 'y must not be empty!');
    nC = size(x,1);
    assert(nC == size(y,1), 'x and y must have same number of channels!');
    assert(round(offset)==offset, 'Offset must be an integer!');
    ly = size(y,2);
    lx = size(x,2);
    prepadding  = func(nC, min(lx, max(0, offset)));
    postpadding = func(nC, max(0, min(lx, lx - offset - ly)));
    ystart = max(1, -offset +1);
    ystopp = min(ly, lx-offset);
    x = x - [prepadding y(:, ystart:ystopp) postpadding];