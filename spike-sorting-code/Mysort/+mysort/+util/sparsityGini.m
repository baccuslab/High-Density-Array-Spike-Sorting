
function s = sparsityGini(c)
% computes the sparsity of vector c after the definition given by Gini.
% If c is a matrix, the sparsity will be computed on the columns of c
% M. O. Lorenz, ÿMethods of measuring concentrations of wealth,ÿ J.
% Amer. Stat. Assoc., 1905

if size(c,1) > 1 && size(c,2) > 1
    s = zeros(1, size(c,2));
    for i=1:size(c,2)
        s(i) = mysort.util.sparsityGini(c(:,i));
    end
    return
end

if ~any(c)
    s = 0;
    return
end

c = sort(c);

N = length(c);

s = 0;
sumc = sum(abs(c));

for k=1:N
    s = s + (c(k)/sumc) * (N-k+.5)/N;
end
s = 1 - 2*s;

