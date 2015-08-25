% function sortNTKflist(flist)
%ana.sorter.initMatlabpool(12);
t_glob = tic();
pd = pdefs();

expname = 'ND10blue';
dpath = fullfile(pd.serverData, 'Antonia', 'config 9', expname);
outpath = fullfile(dpath, 'ntkflistSorter');
if ~exist(outpath, 'file')
    mkdir(outpath);
end

flist = {'wash/Trace_id1140_2013-06-16T15_26_21_14.stream.ntk'
         'wash/Trace_id1140_2013-06-16T15_26_21_15.stream.ntk'
         'wash/Trace_id1140_2013-06-16T15_26_21_16.stream.ntk'
         'wash/Trace_id1140_2013-06-16T15_26_21_17.stream.ntk'
         'drug/Trace_id1140_2013-06-16T15_26_21_18.stream.ntk'
         'drug/Trace_id1140_2013-06-16T15_26_21_19.stream.ntk'
         'drug/Trace_id1140_2013-06-16T15_26_21_20.stream.ntk'
         'drug/Trace_id1140_2013-06-16T15_26_21_21.stream.ntk'
         'washout/Trace_id1140_2013-06-16T15_26_21_22.stream.ntk'
         'washout/Trace_id1140_2013-06-16T15_26_21_23.stream.ntk'
         'washout/Trace_id1140_2013-06-16T15_26_21_24.stream.ntk'
         'washout/Trace_id1140_2013-06-16T15_26_21_25.stream.ntk'};



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
toc(t_glob);
