%% INIT FELIX
% matlabpool(12);
pd = pdefs();
% make sure you have the path of a folder with all .ntk files WITH A SINGLE
% COMMON configuration in the variable "dpath"

depath = fullfile(pd.serverData, 'Katja');

expNames = {'20131210c_ntk'
            '20131217f_ntk'};


        
for kk=1:length(expNames)
    expName = expNames{kk};
    sortingName = 'Sorting1';
    outPath = fullfile(depath, expName, 'Out');
    if ~exist(outPath, 'file'); mkdir(outPath); end
    confListFile = fullfile(depath, expName, 'configNTKLists.mat');
    c = load(confListFile);
    configNTKLists = c.configNTKLists;

    %% Now we CONVERT all ntk files to H5. Doing this we also prefilter the data
    h5FileLists = {};
    for i=1:length(configNTKLists)
        for k=1:length(configNTKLists{i})
            fn = configNTKLists{i}{k};
            [ffolder ffname fsuffix] = fileparts(fn);
            % This is for Katja
            fn = fullfile(depath, expName, fn);
            % end
            fnpf = fullfile(outPath, [ffname fsuffix '_prefilter.h5']);
            h5FileLists{i}{k} = fnpf;
            bConvertFile = true;
            if exist(fnpf, 'file')
                try
                    bIsInProcess = hdf5read(fnpf, '/bFileIsInProcess');
                    if bIsInProcess
                        % This means the file is still in process. We assume no
                        % other user/matlab instance is writing so it is
                        % probably a left over from some other run and we
                        % delete it.
                        delete(fnpf); 
                    else
                        % The file was already properly processed, we dont need
                        % to do anything.
                        bConvertFile = false;
                    end
                catch
                    delete(fnpf);
                end
            end
            if bConvertFile
                % Only convert file if it was non existent or deleted
                mysort.mea.convertNTK2HDF(fn, 'prefilter', 1, 'outFile', fnpf);
            end
        end
    end

    %%
    if 0
        for i=1:length(h5FileLists)
            sortOutPath = fullfile(outPath, ['Config' num2str(i)], sortingName);
%             ana.michele.spikeSortingStartScripts.run(h5FileLists{i}, sortOutPath, sortingName); 
            % SORT
            [gdf_merged T_merged localSorting localSortingID sessionLengths] = ana.startHDSorting(h5FileLists{i}, sortOutPath);
            save(fullfile(sortOutPath, [sortingName '_results.mat']), 'gdf_merged', 'T_merged', 'localSorting', 'localSortingID', 'sessionLengths');
        end
    end

    %% EXPORT TO SINGLE FOLDER
    exportPath = fullfile(outPath, 'Results');
    if ~exist(exportPath, 'file')
        mkdir(exportPath)
    end
    R = {};
    clear RR G
    for i=1:length(h5FileLists)
        sortOutPath = fullfile(outPath, ['Config' num2str(i)], sortingName);
        sourceFile = fullfile(sortOutPath, [sortingName '_results.mat'])
        groupFile = fullfile(fullfile(sortOutPath, 'groupFile.mat'));
        assert(exist(sourceFile, 'file')>0, 'File not found!');
        D = load(sourceFile);
        G(i) = load(groupFile);
        RR(i) = D;

        L = D.sessionLengths;
        % break gdf apart
        start_end_times = [0 cumsum(L)];
        mgdf = mysort.spiketrain.gdf2multiSessionGdf(D.gdf_merged, start_end_times);
        for k=1:length(L)
            gdf = mgdf(mgdf(:,3)==k,[1:2]);
            R{k, i} = gdf;
        end
    end
    save(fullfile(exportPath, [expName '_resultsForKatja.mat']), 'R', 'RR', 'G');
end

%% COPY TO MAIN FOLDER
for kk=1:length(expNames)
    expName = expNames{kk};
    outPath = fullfile(depath, expName, 'Out');
    exportPath = fullfile(outPath, 'Results');
    
    
    configNTKLists = c.configNTKLists;
    inFile = fullfile(exportPath, [expName '_resultsForKatja.mat']);
    outFile = fullfile(depath, [expName '_resultsForKatja.mat']);
    eval(['!cp "' inFile '" "' outFile '"'])
end


