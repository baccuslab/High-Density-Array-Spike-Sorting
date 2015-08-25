%% Store a string in an H5 file
stringToStore = '1234567890';

L = length(stringToStore);
fname = 'test.h5';
h5path = '/testString';
dims = [1 10];
maxDims = [1 10];
h5type = 'H5T_C_S1';
chunkDims = [1 10];
deflation = 0;
bReadOnly = 0;
% Create file and variable
M = mysort.h5.createVariableAndOrFile(fname, h5path, dims, maxDims, h5type, chunkDims, deflation, bReadOnly);
% Store string
M(1,:) = stringToStore;
% Reload String
storedString = M(1,:)

%%
clear M
delete(fname)

