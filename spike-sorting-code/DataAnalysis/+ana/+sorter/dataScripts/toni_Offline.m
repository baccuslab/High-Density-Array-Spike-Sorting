%% INIT FELIX
pd = pdefs();
% make sure you have the path of a folder with all .ntk files WITH A SINGLE
% COMMON configuration in the variable "confPath"
% expname = 'Cortex mice';
expname = 'NHP150416';
sortingName = 'Sorting1';
expPath    = fullfile(pd.serverData, 'Antonia', expname);
outPath    = fullfile(pd.serverData, 'Antonia', [expname '_Out']);
flistFolder = fullfile(expPath, []);
assert(exist(flistFolder, 'file')>0, 'Config Directory not found!')
confLists = dir(fullfile(flistFolder, 'flist*.mat'));

configFolderPathList = {};
configNTKLists = {};
%matlabpool(12)

%% PREPARE ALL CONFIGS AND FILE LISTS
% confNames = {};
fileNotFound = {};
for confi=1:length(confLists)
    fl = load(fullfile(flistFolder, confLists(confi).name));
    flist = fl.flist;
    flistOut = fl.flist;
    
    [bla confName bar] = fileparts(confLists(confi).name);
    confName = confName(6:end); % cut away the useless part of the filename
    confNames{confi} = confName; % important for FRN export
% %     outpath = fullfile(outPath, [confName '_Sorted']);
    confPath = fullfile(expPath, confName);
    confPathOut = fullfile(outPath, confName);
    
    for i=1:length(flist)
        % use expPath here if confPath is already in the flist from Antonia
        flistOut{i} = fullfile(outPath, [flist{i} '_prefilter.h5']);
        flist{i} = fullfile(expPath, flist{i});        
        [outFolder b c] = fileparts(flistOut{i});
        if ~exist(outFolder, 'file')>0
            mkdir(outFolder);
        end
        if ~(exist(flist{i}, 'file')>0)
            fileNotFound{end+1} = flist{i};
        end
    end
    
    configNTKLists{confi} = flist;
    configH5Lists{confi} = flistOut;
    configFolderPathList{confi} = confPathOut;
%     if ~exist(outpath, 'file')
%         mkdir(outpath);
%     end 
end
if ~isempty(fileNotFound)
    logFile = fullfile(pd.serverData, 'Antonia', [expname '_filesNotFound.txt']);
    str = sprintf('File not found! %s\n', fileNotFound{:});
    mysort.util.logToFile(logFile, str);
    assert(isempty(fileNotFound), str)
end

%% Now we CONVERT all ntk files to H5. Doing this we also prefilter the data
if 1
    allFiles = {};
    allFilesOut = {};
    for i=1:length(configNTKLists)
        for j=1:length(configNTKLists{i})
            allFiles{end+1} = configNTKLists{i}{j};
            allFilesOut{end+1} = configH5Lists{i}{j};
        end
    end
    
    parfor combi=1:length(allFiles)
        fn = allFiles{combi};
        fnpf = allFilesOut{combi};
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
%% EXECUTE SORTING
if 1
for i=1:length(configH5Lists)
    logfile = fullfile(configFolderPathList{i}, ['errlog' num2str(i) '.txt']);
    doneFile = fullfile(configFolderPathList{i}, [sortingName '_done'  num2str(i) '.txt']);
    if exist(doneFile, 'file')
        fprintf('This sorting seems to be done. (file exists: %s)\n', doneFile);
        continue
    end
    try
        ana.toni.executeSortingH5(sortingName, configFolderPathList{i}, configH5Lists{i});
        mysort.util.logToFile(doneFile, 'Done.');
    catch
        fprintf('Error sorting config %d!\n', i);
        mysort.util.logLastErrToFile(logfile);
    end    
end
end
resultOutPath = fullfile(outPath, 'Results');

%% EXPORT FrameNumbers for Antonia
if 1
    disp('Starting FRN export');
    frnPath = fullfile(outPath, 'FRN');
    if ~exist(frnPath, 'file')
        mkdir(frnPath);
    end
    for i=1:length(configH5Lists)
        fprintf('config: %d\n', i);
        for k=1:length(configH5Lists{i})
            fn = configNTKLists{i}{k};
            fnpf = configH5Lists{i}{k};
            bIsInProcess = hdf5read(fnpf, '/bFileIsInProcess');
            M = mysort.h5.matrix(fnpf, '/Sessions/Session0/sig');
            % check if binary file
            if min(size(M))==1
                binDims_ = mysort.h5.matrix(fnpf, '/Sessions/Session0/bin_dims');
                binDims = binDims_(:,:);
                [a,b,c] = fileparts(fnpf);
                binFile = fullfile(a, M(:,:));
                M = mysort.ds.binaryFileMatrix(binFile, binDims);
            end
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
ana.toni.exportSorting(sortingName, configFolderPathList, confNames, resultOutPath);
end
%% COPY to NetworkShare
if 1
    sharePath = fullfile(pd.networkTempShare, 'Antonia', expname);
    if ~exist(sharePath, 'file')
        mkdir(sharePath);
    end
    eval(['!cp -R "' resultOutPath '" "' sharePath '"'])
end

%%


[SSC S compoundMea G] = ana.toni.inspectSorting(sortingName, configFolderPathList{maxNeuronConfig}, configNTKLists{maxNeuronConfig});
unitsFound = SSC.unitNames;
nElectrodes = size(compoundMea, 2);


showUnitsIdx = find(unitsFound == maxNeuronID);
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



