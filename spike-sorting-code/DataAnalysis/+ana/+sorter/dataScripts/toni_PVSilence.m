%% INIT FELIX
pd = pdefs();
% make sure you have the path of a folder with all .ntk files WITH A SINGLE
% COMMON configuration in the variable "confPath"
expname = 'PV Silence';
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
            disp(['Converting ' fn]);
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
    if exist(doneFile, 'file')
        continue
    end
    try
        ana.toni.executeSorting(sortingName, configFolderPathList{i}, configNTKLists{i});
        mysort.util.logToFile(doneFile, 'Done.');
    catch
        disp(['Error sorting config %d!' num2str(i)]);
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
    
    % compute modulation index
    SC_ = SC+1;
    Activity         = sum(SC([1:3],:));
    ModulationIndex1 = SC_(2,:)./mean(SC_([1 3],:));
    ModulationIndex2 = SC_(2,:)./SC_(1,:);

    g = struct();
    g.elNr = elNr;
    g.elPos = compoundMea.MultiElectrode.electrodePositions;
    g.elLabels = compoundMea.MultiElectrode.electrodeLabels;
    
    g.goodElectrodes = SI;
    g.goodElectrodesAmplitudes = SM;
    g.mgdf = mgdf;
    g.SingleFileGDF = R;
    g.FileLength = L;
    g.Units = units;
    g.SpikeCountsPerFile = SC;
    g.ModulationIndex1 = ModulationIndex1;
    g.ModulationIndex2 = ModulationIndex2;
    g.TotalSpikeCount  = Activity;
    
    GE(configIdx) = g;
end
toc

%% Compile All Neurons from all Configurations into one big matrix
BM = [];
nEl = 0;
EL = [];
ELNR = [];
for i = 1:length(GE)
    g = GE(i);
    this = [i*ones(length(g.Units),1) g.Units(:) g.TotalSpikeCount(:) g.SpikeCountsPerFile' g.ModulationIndex1(:) g.ModulationIndex2(:) g.goodElectrodes(:,1:20) + nEl g.goodElectrodesAmplitudes(:,1:20)];
    BM = [BM; this];
    EL = [EL; g.elPos];
    ELNR = [ELNR; g.elNr];
    nEl = nEl + size(g.elPos,1);
end
% Columns in sBM represent:
%    1       2         3               4             5            6          7        8             9:28                    29:48
% Config, UnitID, Total#Spikes, Spikes Control, Spikes PSEM, Spikes Wash, ModIdx1, ModIdx2, [20 best Electrodes], [Amplitudes on 20 best ELs]

% Sort BM according to Mod Idx!
sBM = sortrows(BM, -7);
% Remove Neurons with too few spikes
sBM(sBM(:,4)<4,:) = [];
sBM(sBM(:,6)<4,:) = [];
% Remove Neurons with too low amplitude
sBM(sBM(:,30)<20,:) = [];
% Remove Neurons from config 3
sBM(sBM(:,1)==3,:) = [];
% Remove Neurons with much fewer spikes in control than wash
sBM(sBM(:,4)./sBM(:,3) < .33,:) = [];



elIdx = 9:18;
selNeurons = sBM([1:5 end-34:end],:);
ELs4Michele = [];
ELs4Michele(:,:,1) = reshape(EL(selNeurons(:,elIdx),1), [size(selNeurons,1) length(elIdx)]);
ELs4Michele(:,:,2) = reshape(EL(selNeurons(:,elIdx),2), [size(selNeurons,1) length(elIdx)]);
ELs4Michele(:,:,3) = reshape(ELNR(selNeurons(:,elIdx),1), [size(selNeurons,1) length(elIdx)]);
save(fullfile(outPath, 'Electrodes4Michele.mat'), 'ELs4Michele');


%%
% To get config idx and unit name, look into sBM and select a neuron
configIdx = 6;
g = GE(configIdx);
[SSC S compoundMea G] = ana.toni.inspectSorting(sortingName, configFolderPathList{configIdx}, configNTKLists{configIdx});
unitsFound = SSC.unitNames;
cutLeft = 10;
cutLength = 45; 
nElectrodes = size(compoundMea, 2);

%%
showUnitsNames = [1207];
for i=1:length(showUnitsNames)
    
    uIdx = find(unitsFound == showUnitsNames(i));
    uID = unitsFound(uIdx);
    uElGroup = (uID - mod(uID,100))/100;
    spikes = SSC.getWaveforms4UnitIdx(uIdx, cutLeft, cutLength);
    tspikes = mysort.wf.v2t(spikes, nElectrodes);
    EP = compoundMea.MultiElectrode.electrodePositions;
    groupElIdx = G.groupsidx{uElGroup};
    
    mysort.plot.figure([1200 800])
    ah = subplot(1,2,1);
    mysort.plot.waveforms2D(tspikes, EP, 'plotMedian', 1, 'AxesHandle', ah);
    hold on
    plot(EP(groupElIdx,1), EP(groupElIdx,2), 'rx', 'markersize', 14, 'linewidth', 3)
    
    
    ah(2) = subplot(1,2,2);
    hold on
    mm = [inf -inf];
    for f=1:3
        st = g.SingleFileGDF{f};
        st = st(st(:,1)==uID,2);
        if isempty(st)
            continue
        end
        if f==2
            c = 'r';
        else
            c = 'k';
        end    
        plot(st/20000, ones(1, length(st))*f, 'o', 'color', c);
        mm = [min(mm(1), min(st)) max(mm(2), max(st))];
    end
    modulation = g.ModulationIndex1(uIdx);
    change = -(100 - 100*g.ModulationIndex1(uIdx));
    title(['Spike Trains, Modulation: ' sprintf('%.2f', modulation) ' Change: ' sprintf('%.1f',change) '%']);
    set(ah(2), 'ylim', [0 4], 'xlim', [0 mm(2)*1.1]/20000);
    xlabel('Time in seconds');
    ylabel('Condition');
    mysort.plot.figureTitle(['Config: ' num2str(configIdx) ' UnitIdx: ' num2str(uIdx) ' UnitName:' num2str(unitsFound(uIdx))]);
    
end



