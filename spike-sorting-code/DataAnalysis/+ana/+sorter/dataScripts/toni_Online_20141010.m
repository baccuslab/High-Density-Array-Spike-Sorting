%% INIT FELIX
pd = pdefs();
% make sure you have the path of a folder with all .ntk files WITH A SINGLE
% COMMON configuration in the variable "confPath"
% expname = 'Cortex mice';
expname = '141010online';
sortingName = 'Sorting1';
% expPath    = fullfile(pd.mea1kData, 'Antonia', expname);
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
if 1
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
    logfile = fullfile(configFolderPathList{i}, ['errlog' num2str(i) '.txt']);
    doneFile = fullfile(configFolderPathList{i}, ['done'  num2str(i) '.txt']);
%     if exist(doneFile, 'file')
%         fprintf('This sorting seems to be done. (file exists: %s)\n', doneFile);
%         continue
%     end
%     try
        ana.toni.executeSorting(sortingName, configFolderPathList{i}, configNTKLists{i});
        mysort.util.logToFile(doneFile, 'Done.');
%     catch
%         fprintf('Error sorting config %d!', i);
%         mysort.util.logLastErrToFile(logfile);
%     end    
end
end
outPath = fullfile(expPath, 'Results');

%% EXPORT FrameNumbers for Antonia
bDebug = false;
if ~bDebug
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
if ~bDebug
    ana.toni.exportSorting(sortingName, configFolderPathList, confNames, outPath);
end
%% COPY to NetworkShare
if ~bDebug
    sharePath = fullfile(pd.networkTempShare, 'Antonia', expname);
    if ~exist(sharePath, 'file')
        mkdir(sharePath);
    end
    eval(['!cp -R "' outPath '" "' sharePath '"'])
end

%% Extract Good Electrodes for all Templates
clear GE;
disp('Computing good Electrodes For Each Template And Config...');
tic
for configIdx = 1:length(configFolderPathList)
    [SSC S compoundMea G] = ana.toni.inspectSorting(sortingName, configFolderPathList{configIdx}, configNTKLists{configIdx});
    unitsFound = SSC.unitNames;
    cutLeft = 10;
    cutLength = 45; 
    nElectrodes = size(compoundMea, 2);

    % Get spike trains and chop into different files/conditions
    L = compoundMea.X.getAllSessionsLength();
    
    % break gdf apart
    start_end_times = [0 cumsum(L)];
    mgdf = mysort.spiketrain.gdf2multiSessionGdf(S.gdf_merged, start_end_times);
    R = {};
    SC = [];
    units = unique(mgdf(:,1));
    for k=1:length(L)
        gdf = mgdf(mgdf(:,3)==k,[1:2]);
        for u=1:length(units)
            SC(k,u) = length(find(gdf(:,1)==units(u)));
        end
        R{k} = gdf;
    end
    
    % Get Electrode Numbers and best Electrodes for Templates 
    elNr = compoundMea.MultiElectrode.electrodeNumbers;
    T = S.T_merged;
    [mi ma] = mysort.wf.tMinMaxPerTemplate(T);
    mi = abs(mi); ma = abs(ma);
    M = max(mi, ma);
    SI = zeros(size(M));
    SM = SI;
    for i=1:size(SI,1)
        [SM(i,:) SI(i,:)] = sort(M(i,:), 'descend');
    end
    % now we have for each template the electrodes sorted accrording to abs
    % amplitude
    
%     % compute modulation index
    SC_ = SC+1;
    Activity         = sum(SC([1:3],:));
    ModulationIndex1 = SC_(2,:)./mean(SC_([1 3],:));
    ModulationIndex2 = SC_(2,:)./SC_(1,:);
    ModulationIndex3 = SC_(1,:)./SC_(3,:);

    g = struct();
    g.elNr = elNr;
    g.unitsFound = unitsFound;
    g.elPos = compoundMea.MultiElectrode.electrodePositions;
    g.elLabels = compoundMea.MultiElectrode.electrodeLabels;
    
    g.goodElectrodes = SI;
    g.goodElectrodesAmplitudes = SM;
    g.goodElectrodesAmplitudesMax = max(abs(g.goodElectrodesAmplitudes), [], 2);
    g.mgdf = mgdf;
    g.SingleFileGDF = R;
    g.FileLength = L;
    g.Units = units;
    g.SpikeCountsPerFile = SC;
    g.ModulationIndex1 = ModulationIndex1;
    g.ModulationIndex2 = ModulationIndex2;
    g.ModulationIndex3 = ModulationIndex3;
    g.TotalSpikeCount  = Activity;
    
    GE(configIdx) = g;
end
toc

%% Compile All Neurons from all Configurations into one big matrix
BM = [];
nEl = 0;
EL = [];
ELNR = [];
SCall = [];
for i = 1:length(GE)
    g = GE(i);
    %                     1                   2                  3                     4                             5                     6                    7                    8:27                                   28:37
    this = [ones(length(g.Units(:)), 1)*i g.Units(:) g.TotalSpikeCount(:) g.goodElectrodesAmplitudesMax(:) g.ModulationIndex1(:) g.ModulationIndex2(:) g.ModulationIndex3(:) g.goodElectrodes(:,1:20) + nEl g.goodElectrodesAmplitudes(:,1:20)];
    BM = [BM; this];
    EL = [EL; g.elPos];
    nEl = nEl + size(g.elPos,1);
    SCall = [SCall; g.SpikeCountsPerFile'];
    ELNR = [ELNR; g.elNr];
    
end


%% MAKE SUBSELECTION OF NEURONS ACCORDING TO MODULATION INDEX
idxAmp = BM(:,4) > 25; % Amplitude
idxMinSpikeCount = SCall(:,1) > 20;  % Spike Count
idxStable = SCall(:,1) > SCall(:,2) & SCall(:,3) > SCall(:,2); %BM(:,7) > .5 & BM(:,7) < 1.5; % Mod 3
idxMod1up   = BM(:,5) > 1 + .08; % Mod 1
idxMod1down = BM(:,5) < 1 - .08; % Mod 1
FINAL_SELECTION_INDEX = find(idxAmp & idxMod1down & idxStable & idxMinSpikeCount);

nNeuronMax = 30;
[sBM, subIdx] = sortrows(BM(FINAL_SELECTION_INDEX, :), -5);
nBlaIdx = 1:min(nNeuronMax, size(sBM,1));
FINAL_SELECTION_INDEX = FINAL_SELECTION_INDEX(subIdx(nBlaIdx));
sBM = sBM(nBlaIdx,:);


%% PLOT SUB SELECTED NEURONS
mysort.plot.figure([1600 600])
nR = 1;
plots = {'Amp' idxAmp
         'Amp + Stabil' idxAmp & idxStable
         'Amp + Mod1 hoch' idxAmp & idxMod1up
         'Amp + Mod1 runter' idxAmp & idxMod1down
         'Amp + Mod1 runter + Stabil' FINAL_SELECTION_INDEX};
nC = size(plots,1)+1;
ah = zeros(nC, 1);
for i=1:nC-1
    ah(i) = subplot(nR, nC, i);
    idx = plots{i, 2};
    SC_s = SCall(idx, :);
    for k=1:size(SC_s,1)
        hold on
        plot([1 2 3], SC_s(k,:)+1, '.-', 'color', [.5 .5 .5])
    end
    plot([1 2 3], mean(SC_s,1)+1, 'k.-', 'linewidth', 2)
    title(plots{i, 1})
    if i==1
        ylabel('Spike Count');
    end
    hl = legend(sprintf('N = %d', size(SC_s,1)));
    set(hl, 'box', 'off')
    set(ah(i), 'YScale', 'log', 'xlim', [0.5 3.5]);    
end
linkaxes(ah, 'y');

subplot(nR, nC, nR*nC);
SCselected = SCall(plots{end,2}, :);
hist(SCselected(:,2)./mean(SCselected(:,[1 3]), 2), [0:.05:1.2]);
set(gca, 'xlim', [0 1.2])


%%
[maxAmp maxNeuronIdx] = max(BM(:,4));
maxNeuronConfig = BM(maxNeuronIdx, 1);
maxNeuronID = BM(maxNeuronIdx, 2);

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

%%
save(fullfile(outPath, 'Sorting1.mat'), 'sBM', 'BM', 'SCall', 'EL', 'ELNR');

%% EXPORT ELECTRODES OF 40 SELECTED NEURONS FOR CONFIG GENERATION SCRIPT
elIdx = 8:17;
selNeurons = sBM;
ELs4Michele = [];
ELs4Michele(:,:,1) = reshape(EL(selNeurons(:,elIdx),1), [size(selNeurons,1) length(elIdx)]);
ELs4Michele(:,:,2) = reshape(EL(selNeurons(:,elIdx),2), [size(selNeurons,1) length(elIdx)]);
ELs4Michele(:,:,3) = reshape(ELNR(selNeurons(:,elIdx),1), [size(selNeurons,1) length(elIdx)]);
save(fullfile(outPath, 'Electrodes4Michele.mat'), 'ELs4Michele');

%% ROUTE CONFIG
mysort.mea.routeExperimentConfig(outPath, ELs4Michele);
