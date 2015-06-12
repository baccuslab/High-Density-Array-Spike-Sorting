
function c = conv(a, b, shape)
% Compute convolution between a and b with conv2 instead of filter. This is
% the way Matlab2010 is doing it.
if nargin < 3
    shape = 'full';
end
assert(size(a,1) == 1, 'a and b must be row vectors!');
assert(size(b,1) == 1, 'a and b must be row vectors!');
c = conv2(a,b,shape);

if strcmp(shape,'full')
    if length(a) > length(b)
        if size(a,1) == 1 %row vector
            c = c.';
        end
    else
        if size(b,1) == 1 %row vector
            c = c.';
        end
    end
elseif strcmp(shape, 'same')==1 || strcmp(shape, 'valid')==1
    if size(a,1) == 1 %row vector
        c = c.';
    end
else
    error('shape parameter must be full, same or valid!');
end