% function sortNTKflist(flist)
%ana.sorter.initMatlabpool(12);
t_glob = tic();
pd = pdefs();

expname = '20130613_HCact';     
for configIdx = 3:3
    %% 1 AND 2 ARE THROUGH, 3 A FILE WAS MISSING: 'config3/drug/Trace_id1140_2013-06-14T16_27_23_0.stream.ntk'
    
    dpath = fullfile(pd.serverData, 'Antonia', expname);
    
% ['config' num2str(configIdx)]
    confName = sprintf('flist130613Conf_%d', configIdx);
    fl = load(fullfile(dpath, [confName '.mat']));
    flist = fl.(confName);
    outpath = fullfile(dpath, [confName '_Sorted']);
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
    save(fullfile(outpath, [expname '_results.mat']), 'gdfs', 'flist', 'electrodeGroups');
    %%
    disp('Done')
    toc(t_glob);
end
