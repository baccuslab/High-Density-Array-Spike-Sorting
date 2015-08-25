%% INIT FELIX
pd = pdefs();
% make sure you have the path of a folder with all .ntk files WITH A SINGLE
% COMMON configuration in the variable "dpath"
% expName = 'configs_23June2013';
expName = 'configs_2Nov2013';
sortingName = 'Sorting1';
depath = fullfile(pd.networkTempShare, 'Hillier_2013');
outPath = fullfile(depath, [expName 'Out']);
if ~exist(outPath, 'file'); mkdir(outPath); end
confListFile = fullfile(depath, [expName '.mat']);
c = load(confListFile);
configNTKLists = c.configNTKLists;

%% Now we CONVERT all ntk files to H5. Doing this we also prefilter the data
h5FileLists = {};
parfor i=1:length(configNTKLists)
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
exportPath = fullfile(outPath, 'Results');
if ~exist(exportPath, 'file')
    mkdir(exportPath)
end
R = {};
for i=1:length(h5FileLists)
    sortOutPath = fullfile(outPath, ['Config' num2str(i)], sortingName);
    sourceFile = fullfile(sortOutPath, [sortingName '_results.mat'])
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
% 
% 
% %% 
% configIdx = 5;
% [SSC S compoundMea G] = ana.toni.inspectSorting(sortingName, configFolderPathList{configIdx}, configNTKLists{configIdx});
% unitsFound = SSC.unitNames;
% cutLeft = 10;
% cutLength = 45; 
% nElectrodes = size(compoundMea, 2);
% %%
% showUnitsIdx = [2 12 20];
% for i=1:length(showUnitsIdx)
%     uIdx = showUnitsIdx(i);
%     uID = unitsFound(uIdx);
%     uElGroup = (uID - mod(uID,100))/100;
%     spikes = SSC.getWaveforms4UnitIdx(uIdx, cutLeft, cutLength);
%     tspikes = mysort.wf.v2t(spikes, nElectrodes);
%     EP = compoundMea.MultiElectrode.electrodePositions;
%     groupElIdx = G.groupsidx{uElGroup};
%     figure
%     mysort.plot.waveforms2D(tspikes, EP, 'plotMedian', 1);
%     hold on
%     plot(EP(groupElIdx,1), EP(groupElIdx,2), 'rx', 'markersize', 14, 'linewidth', 3)
%     mysort.plot.figureTitle(['UnitIdx: ' num2str(uIdx) ' UnitName:' num2str(unitsFound(uIdx))]);
% end
% 


