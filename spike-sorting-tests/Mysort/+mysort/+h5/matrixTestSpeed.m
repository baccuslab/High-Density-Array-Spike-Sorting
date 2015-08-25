clear M
clear mysort.h5.matrix
% clear all
close all

tmp_path = 'C:\LocalData\Temp\H5MatrixTest\';
fname = 'matrixTest.h5';
fnamemat = [fname '.mat'];

ffile = fullfile(tmp_path, fname);
ffilemat = fullfile(tmp_path, fnamemat);


dims = [100000 100];
nRepetitions = 10;
chunkDims = [100 100];
deflation = 1;
h5type = 'H5T_NATIVE_DOUBLE';

RES_notchunked = [];
for r=1:nRepetitions
    RES_notchunked(r,:) = mysort.h5.matrixTestSpeedHelper(ffile, dims, h5type, [], []);
end


RES_chunked = [];
for r=1:nRepetitions
    RES_chunked(r,:) = mysort.h5.matrixTestSpeedHelper(ffile, dims, h5type, chunkDims, deflation);
end

%%
mysort.plot.figure('w', 1200, 'h', 500);
p=1;
ah(p) = subplot(1,2,p); p=p+1;
names = {'H5create', 'H5write', 'H5read', 'H5readPart', 'matSave', 'matLoad'};
boxplot(RES_notchunked, names)
ylabel('time [s]');
title('not chunked');

ah(p) = subplot(1,2,p); p=p+1;
names = {'H5create', 'H5write', 'H5read', 'H5readPart', 'matSave', 'matLoad'};
boxplot(RES_chunked, names)
title('chunked');

linkaxes(ah, 'y');