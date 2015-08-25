%% INIT FELIX
% matlabpool(12);
pd = pdefs();
% make sure you have the path of a folder with all .ntk files WITH A SINGLE
% COMMON configuration in the variable "dpath"
% expName = 'configs_23June2013';
% expNames = {'configs_15Nov2013'
%             'configs_24July2013'
%             'configs_26Nov2013'};
% 'configs_9Nov2012'
%             'configs_12Nov2012'
%             'configs_20Apr2012'
% expNames = {%'configs_30Nov2012'
%             'configs_29Jun2012'
% expNames = {'configs_10Dec2013'};
% expNames = {'configs_16Dec2013'};
% expNames = {'configs_10Jan2014'};
% expNames = {'configs_11Nov2013'};
%  expNames = {'configs_28Jan2014'};
%  expNames = {'configs_30Jan2014'};
% expNames = {'configs_6Feb2014_A'
%             'configs_6Feb2014_B'
%             'configs_6Feb2014_C'};
% expNames = {'configs_17Feb2014'};
% expNames = {'configs_18Feb2014'};
% expNames = {'configs_4Mar2014_A'
%             'configs_4Mar2014_B'
%             'configs_4Mar2014_C'};
% expNames = {'configs_17Mar2014'
%             'configs_25Mar2014_A'
%             'configs_25Mar2014_B'};
% expNames = {'configs_22Apr2014_A'
%             'configs_22Apr2014_B'};
% expNames = {'configs_3May2014_A'
%             'configs_3May2014_B'};
% expNames = {'configs_27Jun2014_A'
%             'configs_27Jun2014_B'
%             'configs_27Jun2014_C'};
% expNames = {'configs_14Jul2014_A'
%             'configs_14Jul2014_B'};        
% expNames = {'configs_14Jul2014_C'}; 
% expNames = {'configs_14Jul2014_D'}; 
% expNames = {'configs_14Jul2014_G'
%     'configs_14Jul2014_E'
%     'configs_14Jul2014_F'
%     }; 
% expNames = {'configs_26Jul2014_B'
%     'configs_26Jul2014_A'
%     }; 
% expNames = {'configs_29Jun2012_B'
%     }; 
% expNames = {'configs_7Aug2014_A'
%     'configs_7Aug2014_B'
%     'configs_7Aug2014_C'
%     }; 
% expNames = {'configs_8Aug2014_ALL'
%     }; 
% expNames = {'configs_7Aug2014_D'
%             'configs_20Sep2013'
%             'configs_11Oct2013'
%             'configs_16Oct2013'
%     }; 
% expNames = {'configs_10Oct2014_A'
%             'config_10Oct2014_B'};         
% expNames = {'configs_10Oct2014_B'}; 
% expNames = {'configs_13Oct2014_A'}; 
% expNames = {'configs_13Oct2014_B'}; 
% expNames = {'configs_13Oct2014_C'}; 
% expNames = {'configs_15Oct2014_A'
%             'configs_15Oct2014_B'}; 
% expNames = {'configs_20Oct2014_A'
%             'configs_20Oct2014_B'}; 
% expNames = {'configs_31Oct2014_A'
%             'configs_31Oct2014_B'}; 
% expNames = {'configs_3Nov2014_A'}; 
% expNames = {'configs_3Nov2014_B'
%             'configs_3Nov2014_C'}; 
% expNames = {'configs_14Nov2014_A'
%             'configs_14Nov2014_B'
%             'configs_14Nov2014_C'};        
% expNames = {'configs_27Nov2014_A'
%             'configs_27Nov2014_B'}; 
% expNames = {'configs_28Nov2014_A'
%             'configs_28Nov2014_B'
%             'configs_28Nov2014_C'}; 
% expNames = {'configs_27Nov2014_C'
%             'configs_27Nov2014_D'
%             'configs_27Nov2014_E'};
% expNames = {'configs_7Dec2014_A'
%             'configs_7Dec2014_B'};
% expNames = {'configs_7Dec2014_C'};
% expNames = {'configs_7Dec2014_D'
%             'configs_7Dec2014_E'};     
% expNames = {'configs_15Dec2014_A'
%             'configs_15Dec2014_B'
%             'configs_15Dec2014_C'
%             'configs_15Dec2014_D'};
% expNames = {'configs_18Dec2014_A'
%             'configs_18Dec2014_B'
%             'configs_18Dec2014_C'
%             'configs_18Dec2014_D'
%             'configs_18Dec2014_E'};
% expNames = {'configs_18Dec2014_F'};
% expNames = {'configs_18Dec2014_G'};
% 'configs_7Jan2015_A'
% expNames = {'configs_7Jan2015_B'
%             'configs_8Jan2015_A'};
% expNames = {'configs_12Jan2015_A'};
% expNames = {'configs_12Jan2015_B'};
% expNames = {'configs_14Jan2015_A'};
% expNames = {'configs_15Jan2015_A'};
% expNames = {'configs_15Jan2015_B'};
% expNames = {'configs_19Jan2015_A'};
% expNames = {'configs_20Jan2015_A'};
% expNames = {'configs_20Jan2015_B'};
% expNames = {'configs_23Jan2015_A_hamster_center'
%             'configs_22Jan2015_A_monkey'
%             'configs_23Jan2015_B_hamster_periphery'};
% expNames = {'configs_24Jan2015_A'};
% expNames = {'configs_25Jan2015_A'
%             'configs_25Jan2015_B'};
% expNames = {'configs_25Jan2015_C'
%             'configs_25Jan2015_D'};
expNames = {'configs_25Jan2015_D'};


for kk=1:length(expNames)
expName = expNames{kk};
sortingName = 'Sorting2';
depath = fullfile(pd.networkTempShare, 'Hillier_2013');
outPath = fullfile(depath, [expName 'Out']);
if ~exist(outPath, 'file'); mkdir(outPath); end
confListFile = fullfile(depath, [expName '.mat']);
c = load(confListFile);
configNTKLists = c.configNTKLists;

%% Now we CONVERT all ntk files to H5. Doing this we also prefilter the data
allFiles = {};
allCombis = [];
h5FileLists = {};
for i=1:length(configNTKLists)
    for k=1:length(configNTKLists{i})
        fn = configNTKLists{i}{k};
        allFiles{end+1} = fn;
        [ffolder ffname fsuffix] = fileparts(fn);
        fnpf = fullfile(outPath, [ffname fsuffix '_prefilter.h5']);
        h5FileLists{i}{k} = fnpf;
    end
end
   

parfor combi=1:length(allFiles)
    fn = allFiles{combi};
    [ffolder ffname fsuffix] = fileparts(fn);
    fnpf = fullfile(outPath, [ffname fsuffix '_prefilter.h5']);
    bConvertFile = true;
    if exist(fnpf, 'file')
        try
            bIsInProcess = hdf5read(fnpf, '/bFileIsInProcess');
%                 bIsInProcess = 1;
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
%         DS = mysort.mea.CMOSMEA(fnpf);
%         x = DS(1:10,:);
    end
end
fn = [];

%%
for i=1:length(h5FileLists)
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
end
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


