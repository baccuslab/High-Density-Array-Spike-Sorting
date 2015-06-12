cd('/net/bs-filesvr01/export/group/hierlemann/Temp/FelixFranke/LocalData/Michele')
ffile = 'gratings.stream.ntk';

outFile1 = 'test_gratings11.stream.h5';
outFile2 = 'test_gratings21.stream.h5';
outFile3 = 'test_gratings31.stream.h5';
outFile4 = 'test_gratings41.stream.h5';

%% Test first, raw data conversion, no filtering, no downsampling, result will be in uint8
delete(outFile1);
mysort.mea.convertNTK2HDF(ffile, 'outFile', outFile1);

%% Test filtering, no downsampling, result will be in single precision
delete(outFile2);
mysort.mea.convertNTK2HDF(ffile, 'outFile', outFile2, 'prefilter', 1);

%% Test ONLY downsampling, result will be in single precision
delete(outFile3);
mysort.mea.convertNTK2HDF(ffile, 'outFile', outFile3, 'prefilter', 0, 'downSample', 10);

%% Test filtering and downsampling, result will be in single precision
delete(outFile4);
mysort.mea.convertNTK2HDF(ffile, 'outFile', outFile4, 'prefilter', 1, 'downSample', 10, 'hpf', 10, 'lpf', 400);

%% Visualize result, make sure the internal online filter is not used
mea1 = mysort.mea.CMOSMEA(outFile1, 'useFilter', 0, 'name', 'Raw');
mea2 = mysort.mea.CMOSMEA(outFile2, 'useFilter', 0, 'name', 'Prefiltered (300Hz-7000Hz)');
mea3 = mysort.mea.CMOSMEA(outFile3, 'useFilter', 0, 'name', 'ONLY Downsampled (10)');
mea4 = mysort.mea.CMOSMEA(outFile4, 'useFilter', 0, 'name', 'Downsampled AND Prefiltered (10Hz-400Hz)');

mysort.plot.SliderDataAxes({mea1, mea2, mea3, mea4}, 'timeIn', 'sec')

% mysort.plot.SliderDataAxes(mea1, 'timeIn', 'sec')
% X = double(mea1(1:end,1:end));


%% Test to write individual channels into one chunk each. Should make it faster to read a single channel
cd('/net/bs-filesvr01/export/group/hierlemann/Temp/FelixFranke/LocalData/Michele')
% cd('Z:\Temp\FelixFranke\LocalData\Michele')
ffile = 'gratings.stream.ntk';

outFile1 = 'test_gratings_chunkAll.stream.h5';
outFile2 = 'test_gratings_chunkIndividual.stream.h5';
delete(outFile1);
delete(outFile2);

tic
mysort.mea.convertNTK2HDF(ffile, 'deflation', 0, 'outFile', outFile1, 'chunkIndividualChannels', 0);
toc
tic
mysort.mea.convertNTK2HDF(ffile, 'deflation', 0, 'outFile', outFile2, 'chunkIndividualChannels', 1);
toc

%%

tic
X = hdf5read(outFile1, '/Sessions/Session0/sig');
a=toc
clear X
tic
X = hdf5read(outFile2, '/Sessions/Session0/sig');
b=toc
clear X

tic
M1 = mysort.mea.CMOSMEA(outFile1);
toc
M2 = mysort.mea.CMOSMEA(outFile2);

tic
x = M1(:,[1:102]);
c = toc

tic
x = M2(:,[1:102]);
d = toc