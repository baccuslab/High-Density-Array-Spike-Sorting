r = [1 0.9 0.8];

C = toeplitz(r);
cond(C)

D = mysort.util.diagonalSubspaceLoading(C, 10);

cond(D)

D = mysort.util.diagonalSubspaceLoading(C, 15.3);

cond(D)