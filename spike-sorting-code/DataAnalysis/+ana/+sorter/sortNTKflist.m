% function sortNTKflist(flist)
ana.sorter.initMatlabpool(12);
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
% flist = {'Trace_id858_2012-04-19T09_49_10_26.stream.ntk'
%          'Trace_id858_2012-04-19T09_49_10_27.stream.ntk'
%          'Trace_id858_2012-04-19T09_49_10_28.stream.ntk'
%          'Trace_id858_2012-04-19T09_49_10_29.stream.ntk'
%          'Trace_id858_2012-04-19T09_49_10_30.stream.ntk'
%          'Trace_id858_2012-04-19T09_49_10_31.stream.ntk'};
flist = flist(10:end);

%% INIT     
if 0
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
