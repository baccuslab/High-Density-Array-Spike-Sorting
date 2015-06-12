
function s = sparsity(x, mode)
% computes the sparsity of vector x

if nargin == 1
    mode = 'Gini';
end

if strcmp(mode, 'Gini');
    s = mysort.util.sparsitiyGini(x);
elseif strcmp(mode, 'Hoyer');
    s = mysort.util.sparsitiyHoyer(x);
end

