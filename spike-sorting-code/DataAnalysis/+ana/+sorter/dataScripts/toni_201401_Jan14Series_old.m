% function sortNTKflist(flist)
% ana.sorter.initMatlabpool(12);
clear all
t_glob = tic();
pd = pdefs();

expname = 'Jan 14 series';
depath    = fullfile(pd.serverData, 'Antonia', expname);
flistpath = fullfile(depath, ['flists ' expname]);
assert(exist(flistpath, 'file')>0, 'Config Directory not found!')
confLists = dir(fullfile(flistpath, 'flist*.mat'));

% Test all config lists first
if 0
    for confi=1:length(confLists)
        fl = load(fullfile(flistpath, confLists(confi).name));
        flist = fl.flist;
        ntksorter = ana.sorter.ntkFileListSorter(depath, flist, '');
        V = ntksorter.checkFList();
    end
end

%%
for confi=1:length(confLists)
    % stopped in 8 because of 1 sample long chunk! [4 7 8 9] 
    % that is resolved now. But problem in 9 with flist:
    % config9/white_noise/control/Trace_id1140_2013-06-16T15_26_21_7.stream.ntk
    
    fl = load(fullfile(flistpath, confLists(confi).name));
    flist = fl.flist;
    [bla confName bar] = fileparts(confLists(confi).name);
    confName = confName(6:end);
    
    outpath = fullfile(depath, [confName '_Sorted']);
    dpath = depath;
    if ~exist(outpath, 'file')
        mkdir(outpath);
    end

    %% INIT     
    if 1
        disp('Start')
        tic
        ntksorter = ana.sorter.ntkFileListSorter(dpath, flist, outpath);
        toc
        V = ntksorter.checkFList();
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
    templates = {};
    for i=1:length(ntksorter.singleElGroupSorters)
        templates{i} = ntksorter.singleElGroupSorters(i).templates;
    end
    save(fullfile(outpath, [expname '_results.mat']), 'gdfs', 'flist', 'electrodeGroups', 'templates');
    %%
    disp('Done')
    toc(t_glob);
end
