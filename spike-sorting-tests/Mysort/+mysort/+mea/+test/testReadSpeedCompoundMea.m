% pd = pdefs;
% % dpath = 'C:\LocalData\Antonia\H5ReadTest';
% dpath = '/net/bs-filesvr01/export/group/hierlemann/recordings/collaborations/Antonia/Nov13 series/131105conf1/ctrl';
% 
% fnames = dir(fullfile(dpath, '*.h5'));
% filth5fnames = {};
% for i=1:length(fnames)
%     filth5fnames{i} = fullfile(dpath, fnames(i).name);
% end
% 
% %%
% compoundMea = mysort.mea.compoundMea(filth5fnames, 'useFilter', 0, 'name', 'PREFILT');
% 
% %%
tic
X = hdf5read(filth5fnames{1}, '/Sessions/Session0/sig');
toc
% % 17.86s on Laptop
% % 5.9 on GPU03
tic
X = double(hdf5read(filth5fnames{1}, '/Sessions/Session0/sig'))*100.1;
toc
% % ?s on Laptop
% % 6.9 on GPU03
%%
nS = size(compoundMea.X.sessionList(1), 1);
% nS = 10;
tic
X = compoundMea(1:nS, :);
toc
% 11.31s on GPU03

%%
times = [100 : 100 :nS];
cutLeft = 20;
cutLength = 50;
tic
wfs = compoundMea.getWaveform(times, cutLeft, cutLength);
toc
