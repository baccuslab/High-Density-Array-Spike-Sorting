
function [i, j, val] = matrixArgMax(M)
    % Computes the argmax in (row, col) of a matrix
    [m, i] = max(M,[],1);
    [val, j] = max(m,[],2);
    i=i(j);