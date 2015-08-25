% function sortNTKflist(flist)
%ana.sorter.initMatlabpool(12);
t_glob = tic();
pd = pdefs();

expname = '20130616_HCinh';     
for configIdx = [9]
    % stopped in 8 because of 1 sample long chunk! [4 7 8 9] 
    % that is resolved now. But problem in 9 with flist:
    % config9/white_noise/control/Trace_id1140_2013-06-16T15_26_21_7.stream.ntk
    
    dpath = fullfile(pd.serverData, 'Antonia', expname);
    
    confName = sprintf('flist130616Conf_%d', configIdx);
    fl = load(fullfile(dpath, [confName '.mat']));
    flist = fl.(confName);
    outpath = fullfile(dpath, [confName '_Sorted']);
    if ~exist(outpath, 'file')
        mkdir(outpath);
    end

    if configIdx == 9
        % make manual fix for flist of config 9
        flist = flist(1:8);
        flist{end+1} = 'config9/white_noise/drug/Trace_id1140_2013-06-16T15_26_21_8.stream.ntk';
        flist{end+1} = 'config9/white_noise/drug/Trace_id1140_2013-06-16T15_26_21_9.stream.ntk';
        flist{end+1} = 'config9/white_noise/wash/Trace_id1140_2013-06-16T15_26_21_7.stream.ntk';
        flist{end+1} = 'config9/white_noise/washout/Trace_id1140_2013-06-16T15_26_21_10.stream.ntk';
        flist{end+1} = 'config9/white_noise/washout/Trace_id1140_2013-06-16T15_26_21_11.stream.ntk';
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
