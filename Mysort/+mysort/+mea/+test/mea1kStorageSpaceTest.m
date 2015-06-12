% This is very weird. It seems that writing UINT16 is faster (>50%) than
% writing INT16 although the file has the same size and the data was in
% INT16... 

X = int16(round(1000*abs(randn(100000, 900))));



%%
testFile = 'testFile1.h5';
deflation = [];
chunkDims = [];
maxDims = size(X);
dims = size(X);
h5Type = 'H5T_NATIVE_USHORT';
if exist(testFile, 'file'); clear sig; delete(testFile); end
tic
sig = mysort.h5.createVariableAndOrFile(testFile, '/Sessions/Session0/sig', ...
    dims, maxDims, h5Type, chunkDims, deflation);
sig(:,:) = X;
clear sig
toc

pause(1)
%%
testFile = 'testFile2.h5';
h5Type = 'H5T_NATIVE_SHORT';
if exist(testFile, 'file'); clear sig; delete(testFile); end
tic
sig = mysort.h5.createVariableAndOrFile(testFile, '/Sessions/Session0/sig', ...
    dims, maxDims, h5Type, chunkDims, deflation);
sig(:,:) = X;
clear sig
toc