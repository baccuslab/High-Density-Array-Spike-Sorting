%% INIT FELIX
% matlabpool(12);
pd = pdefs();
% make sure you have the path of a folder with all .ntk files WITH A SINGLE
% COMMON configuration in the variable "dpath"

% PUT HERE THE !!!NAME OF THE file WITHOUT .mat!!! that contains the config definition
% in a variable called "configNTKLists"
expNames = {'2013_12_11_1281_configs'};
% config configNTKLists is a cell array, each cell represents one config.
% In each of those cells is another cell array with the FULL FILE NAMES of
% all .ntk files (FULL meaning the whole path!)

% THIS IS THE ROOT PATH WHERE ALL THE OUTPUT OF THE SORTINGS AND THE H5 FILES GO
% THIS IS ALSO THE FOLDER WHERE YOU NEED TO PUT THE .mat FILE CONTAINING
% THE CONFIGS (name specified above in expNames).
depath = '/net/bs-filesvr01/export/group/hierlemann/recordings/SpikeSortingOut/Wei/';

% THIS SHOULD START THE SORTING
for kk=1:length(expNames)
    fprintf('Running experiment %s\n', expNames{kk});
    expName = expNames{kk};
    sortingName = 'Sorting1';
    outPath = fullfile(depath, [expName 'Out']);
    if ~exist(outPath, 'file'); mkdir(outPath); end
    confListFile = fullfile(depath, [expName '.mat']);
    assert(exist(confListFile, 'file') > 0, 'Config File not found %s!', confListFile);
    c = load(confListFile);
    configNTKLists = c.configNTKLists;
    nConfigs = length(configNTKLists);
    fprintf('Found Configs: %d\n', nConfigs);
    assert(nConfigs > 0, 'No Config in Config File!');
    %% Now we CONVERT all ntk files to H5. Doing this we also prefilter the data
    h5FileLists = {};
    for i=1:nConfigs
        for k=1:length(configNTKLists{i})
            fn = configNTKLists{i}{k};
            [ffolder ffname fsuffix] = fileparts(fn);
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
    parfor i=1:length(h5FileLists)
        sortOutPath = fullfile(outPath, ['Config' num2str(i)], sortingName);
        ana.michele.spikeSortingStartScripts.run(h5FileLists{i}, sortOutPath, sortingName); 
    end

    %% EXPORT TO SINGLE FOLDER
    R = {};
    for i=1:length(h5FileLists)
        sortOutPath = fullfile(outPath, ['Config' num2str(i)], sortingName);
        sourceFile = fullfile(sortOutPath, [sortingName '_results.mat']);
        assert(exist(sourceFile, 'file')>0, 'File not found!');
        D = load(sourceFile);

        % get individual files' length in samples
        compoundMea = mysort.mea.compoundMea(h5FileLists{i}, 'useFilter', 0, 'name', 'PREFILT');
        L = compoundMea.X.getAllSessionsLength();

        % break gdf apart
        start_end_times = [0 cumsum(L)];
        mgdf = mysort.spiketrain.gdf2multiSessionGdf(D.gdf_merged, start_end_times);
        for k=1:length(L)
            gdf = mgdf(mgdf(:,3)==k,[1:2]);
            R{k, i} = gdf;
        end
    end
    save(fullfile(depath, [expName '_resultsForMichele.mat']), 'R');


    %% EXPORT FOR WEI
    clear neuronInformation;
    for i=1:length(h5FileLists)
        sortOutPath = fullfile(outPath, ['Config' num2str(i)], sortingName);
        sourceFile = fullfile(sortOutPath, [sortingName '_results.mat']);
        groupFile = fullfile(sortOutPath, ['groupFile.mat']);
        assert(exist(sourceFile, 'file')>0, 'File not found!');
        D = load(sourceFile);
        G = load(groupFile);

        ni = struct();
        ni.neuronIDs = unique(D.gdf_merged(:,1));
        nN = length(ni.neuronIDs);

        ni.spikeCounts = zeros(1,nN);
        for k=1:nN
            ni.spikeCounts(k) = sum(D.gdf_merged(:,1)==ni.neuronIDs(k));
        end

        ni.neuronPositions = zeros(2, nN);
        ni.neuronAmplitudes = zeros(1, nN);
        [mi ma mi_idx ma_idx] = mysort.wf.tMinMaxPerTemplate(D.T_merged);
        for k=1:nN
            [mii minElIdx] = min(mi(k,:));
            ni.neuronPositions(:,k) = G.electrodePositions(minElIdx,:);
            ni.neuronAmplitudes(k) = mii;
        end    
        ep = G.electrodePositions;
        ni.electrodePositions = ep;
        ni.electrodeBB_mix_miy_max_may = [min(ep) max(ep)];

        neuronInformation(i) = ni;
        
        if 0
            %%
            mysort.plot.templates2D(D.T_merged, G.electrodePositions, 9)
            plot(ni.neuronPositions(1,:), ni.neuronPositions(2,:), 'kx', 'linewidth', 2, 'markersize', 14)
        end

    end
    save(fullfile(depath, [expName '_resultsForWei.mat']), 'neuronInformation');

end
