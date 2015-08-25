% function sortNTKflist(flist)
%ana.sorter.initMatlabpool(12);
t_glob = tic();
pd = pdefs();
expname = '19Apr2012_DSGCs_rabbit';
cd(fullfile(pd.svnDataRoska, expname, 'Matlab'))
dpath = fullfile(pd.serverDataRoska, expname, 'proc');
outpath = fullfile(pd.sortingOutPath, expname, 'ntkflistSorter');
if ~exist(outpath, 'file')
    mkdir(outpath);
end
flist = {};
flist_rabbit_ds();


%% INIT     
if 1
    disp('Start')
    tic
    ntksorter = ana.sorter.ntkFileListSorter(dpath, flist, outpath);
    toc
    disp('Load Data')
    tic
    ntksorter.loadDataForTraining();
    toc
    disp('Train')
    tic
    ntksorter.train();
    toc
    disp('Save')
    tic
    ntksorter.save();
    toc
else
  %%
    disp('Load Sorter')
    tic
    ntksorter = ana.sorter.ntkFileListSorter(dpath, flist, outpath);
    ntksorter.load();
    toc  
end

%% TM
disp('TM')
tic
ntksorter.templateMatching();
disp('TM')
toc
gdfs = ntksorter.gdfs;
electrodeGroups = ntksorter.electrodes;
save(fullfile(outpath, [expname '_results.mat']), 'gdfs', 'flist', 'electrodeGroups');
%%
disp('Done')
