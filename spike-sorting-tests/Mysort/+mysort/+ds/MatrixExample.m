clear all
close all
clear all
X = randn(100,4);

% create a datasource from a simple matrix
M = mysort.ds.Matrix(X);

size(M)
M(1:10, end-1:end)

M.getData(1:3, 2:3)

M.getSampleInterval()
