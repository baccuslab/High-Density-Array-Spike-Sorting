mysort.util.buildCombinatorics([2 1 2])

% faster:
clear A;
[A{1:3}] = ndgrid([1:2], 1, [1:2]);
A = reshape(cat(3+1,A{:}),[],3)