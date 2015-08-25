%% INIT FELIX
pd = pdefs();
% make sure you have the path of a folder with all .ntk files WITH A SINGLE
% COMMON configuration in the variable "confPath"
expname = 'HCs April 14';
sortingName = 'Sorting1';
expPath    = fullfile(pd.serverData, 'Antonia', expname);
flistFolder = fullfile(expPath, []);
assert(exist(flistFolder, 'file')>0, 'Config Directory not found!')
confLists = dir(fullfile(flistFolder, 'flist*.mat'));

configFolderPathList = {};
configNTKLists = {};
%matlabpool(12)

%% PREPARE ALL CONFIGS AND FILE LISTS
confNames = {};
for confi=1:length(confLists)
    fl = load(fullfile(flistFolder, confLists(confi).name));
    flist = fl.flist;
    
    [bla confName bar] = fileparts(confLists(confi).name);
    confName = confName(6:end); % cut away the useless part of the filename
    confNames{confi} = confName;
    outpath = fullfile(expPath, [confName '_Sorted']);
    confPath = fullfile(expPath, confName);
    
    for i=1:length(flist)
        % use expPath here if conPath is already in the flist from Antonia
        flist{i} = fullfile(expPath, flist{i});
        assert(exist(flist{i}, 'file')>0, sprintf('File not found! %s', flist{i}))
    end
    
    configNTKLists{confi} = flist;
    configFolderPathList{confi} = confPath;
%     if ~exist(outpath, 'file')
%         mkdir(outpath);
%     end 
end   

%% Now we CONVERT all ntk files to H5. Doing this we also prefilter the data
if 0
parfor i=1:length(configNTKLists)
    for k=1:length(configNTKLists{i})
        fn = configNTKLists{i}{k};
        fnpf = [fn '_prefilter.h5'];
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
end
%% EXECUTE SORTING
if 1
for i=1:length(configNTKLists)
    logfile = ['errlog' num2str(i) '.txt'];
    doneFile = ['done'  num2str(i) '.txt'];
    if exist(doneFile, 'file')
        continue
    end
    try
        ana.toni.executeSorting(sortingName, configFolderPathList{i}, configNTKLists{i});
        mysort.util.logToFile(doneFile, 'Done.');
    catch
        disp('Error sorting config %d!', i);
        mysort.util.logLastErrToFile(logfile);
    end    
end
end
outPath = fullfile(expPath, 'Results');

%% EXPORT FrameNumbers for Antonia
if 1
    disp('Starting FRN export');
    frnPath = fullfile(outPath, 'FRN');
    if ~exist(frnPath, 'file')
        mkdir(frnPath);
    end
    for i=1:length(configNTKLists)
        fprintf('config: %d\n', i);
        for k=1:length(configNTKLists{i})
            fn = configNTKLists{i}{k};
            fnpf = [fn '_prefilter.h5'];
            bIsInProcess = hdf5read(fnpf, '/bFileIsInProcess');
            M = mysort.h5.matrix(fnpf, '/Sessions/Session0/sig');
            frn = M(:, 131);
            filename = fn;
            save(fullfile(frnPath, [confNames{i} '_file' num2str(k) '_frn.mat']), 'frn', 'filename');            
        end
    end
    disp('Done.');
end
% return
%% EXPORT TO SINGLE FOLDER
if 1
ana.toni.exportSorting(sortingName, configFolderPathList, confNames, outPath);
end
%% COPY to NetworkShare
if 1
    sharePath = fullfile(pd.networkTempShare, 'Antonia', expname);
    if ~exist(sharePath, 'file')
        mkdir(sharePath);
    end
    eval(['!cp -R "' outPath '" "' sharePath '"'])
end

%% 
configIdx = 7;
[SSC S compoundMea G] = ana.toni.inspectSorting(sortingName, configFolderPathList{configIdx}, configNTKLists{configIdx});
unitsFound = SSC.unitNames;
cutLeft = 10;
cutLength = 45; 
nElectrodes = size(compoundMea, 2);
%%
showUnitsIdx = [2 12 20];
for i=1:length(showUnitsIdx)
    uIdx = showUnitsIdx(i);
    uID = unitsFound(uIdx);
    uElGroup = (uID - mod(uID,100))/100;
    spikes = SSC.getWaveforms4UnitIdx(uIdx, cutLeft, cutLength);
    tspikes = mysort.wf.v2t(spikes, nElectrodes);
    EP = compoundMea.MultiElectrode.electrodePositions;
    groupElIdx = G.groupsidx{uElGroup};
    figure
    mysort.plot.waveforms2D(tspikes, EP, 'plotMedian', 1);
    hold on
    plot(EP(groupElIdx,1), EP(groupElIdx,2), 'rx', 'markersize', 14, 'linewidth', 3)
    mysort.plot.figureTitle(['UnitIdx: ' num2str(uIdx) ' UnitName:' num2str(unitsFound(uIdx))]);
end



