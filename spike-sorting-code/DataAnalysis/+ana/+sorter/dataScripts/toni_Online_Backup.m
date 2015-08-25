%% INIT 
pd = pdefs();
% make sure you have the path of a folder with all .ntk files WITH A SINGLE
% COMMON configuration in the variable "confPath"
% expname = 'Cortex mice';
expname = '150630'; %'150519offline'; %150512wnoi %150512offline % 150519bscreen
sortingName = 'Sorting1';
expDataPath          = fullfile(pd.mea1kData, 'Antonia', expname);
intermediateDataPath = fullfile(pd.mea1kIntermediate, 'Antonia', expname);
resultDataPath       = fullfile(pd.serverData, 'Antonia', expname);

h5Path = [];

% NEW VERSION FOR MEA 1K
flist = {};
h5Files = dir(fullfile(expDataPath, '*.h5'));
nrkFiles = dir(fullfile(expDataPath, 'config', '*.mapping.nrk'));
for i=1:length(h5Files)
    flist{i,1} = fullfile(expDataPath, h5Files(i).name);
    flist{i,2} = fullfile(expDataPath, 'config', nrkFiles(1).name);
end


if ~exist(intermediateDataPath, 'file')
    mkdir(intermediateDataPath);
end
if ~exist(resultDataPath, 'file')
    mkdir(resultDataPath);
end

configFolderOutPathList = {};
configNTKLists = {};
%matlabpool(12)

%% PREPARE ALL CONFIGS AND FILE LISTS
confNames = {};
configH5Lists = {};
for confi=1:1

    confName = expname;
    confNames{confi} = confName;
    confOutPath = fullfile(intermediateDataPath, confName);
    if ~exist(confOutPath, 'file')
        mkdir(confOutPath);
    end
    
    configNTKLists{confi} = flist;
    configFolderOutPathList{confi} = confOutPath;
    configH5Lists{confi} = [];
    for k=1:length(configNTKLists{confi})
        [a, b, c] = fileparts(configNTKLists{confi}{k,1});
        configH5Lists{confi}{k} = fullfile(confOutPath, [b c '_prefilter.h5']);        
    end
%     if ~exist(outpath, 'file')
%         mkdir(outpath);
%     end 
end   

%% Now we CONVERT all ntk files to H5. Doing this we also prefilter the data
if 1
    for i=1:length(configNTKLists)
        parfor k=1:length(configNTKLists{i})
            fn   = configNTKLists{i}{k,1};
            fnpf = configH5Lists{i}{k};
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
                X = dlmread(configNTKLists{i}{k,2});
                isNotConnected = any(X==-1,2);
                X(isNotConnected,:) = [];
                % Im H5 file ist spalte 0 channel 0.
                % 1023 ist der letzte Kanal dann DAC und !command counter"
                % DAC ist nicht im nrk file 
                idx = X(:,1) + 1;
                m = struct();
                m.mposx(idx) = X(:,2);
                m.mposy(idx) = X(:,3);
                m.elNo(idx)   = X(:,1);
                m.chans(idx)  = X(:,1);
                m.noChans = length(m.elNo);
                m.file = fname;         
                
                % Only convert file if it was non existent or deleted
                mysort.mea.preprocessMea1kH5File(fn, m, 'prefilter', 1, 'outFile', fnpf, 'h5path', h5Path);
            end
        end
    end
end
%% EXECUTE SORTING
if 1
    for i=1:length(configNTKLists)
        logfile = fullfile(configFolderOutPathList{i}, ['errlog' num2str(i) '.txt']);
        doneFile = fullfile(configFolderOutPathList{i}, ['done'  num2str(i) '.txt']);
    %     if exist(doneFile, 'file')
    %         fprintf('This sorting seems to be done. (file exists: %s)\n', doneFile);
    %         continue
    %     end
    %     try
            ana.toni.executeSorting(sortingName, configFolderOutPathList{i}, configH5Lists{i});
            mysort.util.logToFile(doneFile, 'Done.');
    %     catch
    %         fprintf('Error sorting config %d!', i);
    %         mysort.util.logLastErrToFile(logfile);
    %     end    
    end
end
resultFolder = fullfile(resultDataPath, 'Results');

%% EXTRACT FRAMENUMBERS
logFolderNames = {'control', 'PSEM', 'wash'};

clear STIMINFO;
for i=1:length(configNTKLists)
    for k=1:length(configNTKLists{i})
        rawh5File = configNTKLists{i}{k};
        h5File  = configH5Lists{i}{k};
        stimFile = [rawh5File(1:end-6) 'mat'];
        % cut away raw.h5_prefilter.h5
        cmdFile = [rawh5File(1:end-2) 'cmd'];

        firstFrameNumber = double(hdf5read(h5File, '/Sessions/Session0/frame_numbers/first_fn'));
        CMD = load(cmdFile);
        StimInfo = load(stimFile);
        
%         CMOS = mysort.mea.CMOSMEA(configH5Lists{i}{k,1});
        clear SI;
        SI.stimInfoIndex = [];
        SI.counter = [];
        SI.cmd_index = [];
        SI.framenumber = [];
        SI.sampleTime = []; 
        SI.endSampleTime = [];
        for j=1:length(StimInfo.stimulus_frame_info)
            if (strcmp(StimInfo.stimulus_frame_info{j}.stimulus_type, 'show_image') || ...
                strcmp(StimInfo.stimulus_frame_info{j}.stimulus_type, 'show_checkerboard')) && ...
               ~StimInfo.stimulus_frame_info{j}.is_last
               counter   = double(StimInfo.stimulus_frame_info{j}.counter);
               cmd_index = counter*2 + 1;
               framenumber = double(CMD(cmd_index, 1));
               sampleTime  = framenumber - firstFrameNumber + 1;
               
               SI.stimInfoIndex = [SI.stimInfoIndex j];
               SI.counter = [SI.counter counter];
               SI.cmd_index = [SI.cmd_index cmd_index];
               SI.framenumber = [SI.framenumber framenumber];
               SI.startSampleTime = [SI.sampleTime sampleTime];
               
               counter   = double(StimInfo.stimulus_frame_info{j+1}.counter);
               cmd_index = counter*2 + 1;
               framenumber = double(CMD(cmd_index, 1));
               sampleTime  = framenumber - firstFrameNumber + 1;
               SI.endSampleTime = [SI.endSampleTime sampleTime];
               
               STIMINFO(k,i) = SI;
            end
        end
    end
end

%% EXPORT TO SINGLE FOLDER
if ~bDebug
    ana.toni.exportSorting(sortingName, configFolderOutPathList, confNames, resultFolder);
end



%% Extract Good Electrodes for all Templates
clear GE;
disp('Computing good Electrodes For Each Template And Config...');
tic
for configIdx = 1:length(configFolderOutPathList)
    [SSC S compoundMea G] = ana.toni.inspectSorting(sortingName, configFolderOutPathList{configIdx}, configH5Lists{configIdx});
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
    R_stimuli = {};
    SC = [];
    units = unique(mgdf(:,1));
    for k=1:length(L)
        gdf = mgdf(mgdf(:,3)==k,[1:2]);
        for u=1:length(units)
            SC(k,u) = length(find(gdf(:,1)==units(u)));
        end
        R{k} = gdf;
        TS = [STIMINFO(k,configIdx).startSampleTime(:) STIMINFO(k,configIdx).endSampleTime(:)];
        R_stimuli{k} = mysort.spiketrain.gdf2multiSessionGdf(gdf, TS, [], 0);
    end
    
    %% TODO: ONLY DO THIS FOR WHITE NOISE
    WHITENOISE_INDEX = 4;
    HIST_LENGTH = 11;
    rawh5File = configNTKLists{configIdx}{WHITENOISE_INDEX};
    stimFile = [rawh5File(1:end-6) 'mat'];
    StimInfo = load(stimFile);    
    
    gdf = R_stimuli{WHITENOISE_INDEX};
    stimStartTime = STIMINFO(WHITENOISE_INDEX, configIdx).startSampleTime;
    WNST = StimInfo.stimulus_frame_info{1, STIMINFO(WHITENOISE_INDEX, configIdx).stimInfoIndex}.parameters.colors;
    
    U = unique(gdf(:,1));
    RF = zeros(size(WNST,2), size(WNST,3), HIST_LENGTH, length(U));
    for ui = 1:length(U)
        st = gdf(gdf(:,1) == U(ui),2);
        % compute spike times as index into movie frame
        framerate = 60; % Hz
        st = round(st/20000 * framerate);
        
        for f = 1:HIST_LENGTH
            RF(:,:,f,ui) = mean(WNST(st-(f-1),:,:));
        end        
    end    
        
    
    %% Get Electrode Numbers and best Electrodes for Templates 
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
    %% TODO: ONLY DO THIS FOR PSEM WASH & CONTROL
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
    
    g.SingleStimulusGDF = R_stimuli;
    
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
idxMinSpikeCount = SCall(:,1) > 80;  % Spike Count
idxStable = SCall(:,1) > SCall(:,2); %& SCall(:,3) > SCall(:,2); %BM(:,7) > .5 & BM(:,7) < 1.5; % Mod 3
idxMod1up   = BM(:,5) > 1 + .08; % Mod 1
idxMod1down = BM(:,5) < 2 - .08; % Mod 1
FINAL_SELECTION_INDEX = find(idxAmp & idxMod1down & idxStable & idxMinSpikeCount);

nNeuronMax = 30;
[sBM, subIdx] = sortrows(BM(FINAL_SELECTION_INDEX, :), 5);
nBlaIdx = 1:min(nNeuronMax, size(sBM,1));
sFINAL_SELECTION_INDEX = FINAL_SELECTION_INDEX(subIdx(nBlaIdx));
sBM = sBM(nBlaIdx,:);


%% PLOT SUB SELECTED NEURONS
mysort.plot.figure([1600 600])
nR = 1;
plots = {'Amp' idxAmp
         'Amp + Stabil' idxAmp & idxStable
         'Amp + Mod1 hoch' idxAmp & idxMod1up
         'Amp + Mod1 runter' idxAmp & idxMod1down
         'Amp + Mod1 runter + Stabil' sFINAL_SELECTION_INDEX};
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

ah(end) = subplot(nR, nC, nR*nC);
SCselected = SCall(plots{end,2}, :);
hist(SCselected(:,2)./mean(SCselected(:,[1 3]), 2), [0:.05:1.2]);
title('Contr. vs Silence Ratio');
ylabel('Count');
set(gca, 'xlim', [0 1.2])

%% PLOT SOME SELECTED NEURONS
window = .15;
sampling = .05;
timePeriod = [0 7];
for i = 1:min(30, length(sFINAL_SELECTION_INDEX));
    config = BM(sFINAL_SELECTION_INDEX(i), 1);
    UID    = BM(sFINAL_SELECTION_INDEX(i), 2);
    
    mysort.plot.figure([1200 900]);
    ah = subplot(5,3,[1 2 4 5]);
    hold on
    ah(2) = subplot(5,3, [7 8 10 11]);
    hold on
    ah(3) = subplot(5,3, [13 14]);
    hold on
    FR = [];
    for k=1:3
        for j=1:5
            offset = 4 - (k + j/7);
            idx = GE(config).SingleStimulusGDF{k}(:,1) == UID & ...
                  GE(config).SingleStimulusGDF{k}(:,3) == j;
            ts =  GE(config).SingleStimulusGDF{k}(idx,2);
            [fr_, x] = mysort.spiketrain.toFiringRate(ts/20000, window, sampling, timePeriod);
            FR(j,:) = fr_(:)';
            plot(ah(1), ts/20000, offset*ones(1,length(ts)), '.', 'color', mysort.plot.vectorColor(k))
        end
        plot(ah(2), x, mean(FR,1), '.-', 'color', mysort.plot.vectorColor(k), 'linewidth', 2)
        plot(ah(3), x, std(FR,[],1)./mean(FR,1), '.-', 'color', mysort.plot.vectorColor(k), 'linewidth', 2)
    end
    
    ah(4) = subplot(5,3,3);
    clear bar
    bar([1 2 3], SCall(sFINAL_SELECTION_INDEX(i),:));
    
    ah(5) = subplot(5,3,6);
%     bar([1 2 3 4 5], 
    
    hl = legend(ah(2), 'control', 'PSEM', 'wash');
    set(hl, 'box', 'off')
    xlabel(ah(3), 'Time [sec]');
    ylabel(ah(2), 'Firing Rate');
    ylabel(ah(1), 'Trials / Conditions');
    str = sprintf('exp: %s, cellidx: %d, cell unit: %d', expname, sFINAL_SELECTION_INDEX(i), sBM(i,2));
    mysort.plot.figureTitle(str)
    linkaxes(ah, 'x')
end
%%
[maxAmp maxNeuronIdx] = max(BM(:,4));
maxNeuronConfig = BM(maxNeuronIdx, 1);
maxNeuronID = BM(maxNeuronIdx, 2);

[SSC S compoundMea G] = ana.toni.inspectSorting(sortingName, configFolderOutPathList{maxNeuronConfig}, configH5Lists{maxNeuronConfig});
unitsFound = SSC.unitNames;
nElectrodes = size(compoundMea, 2);

idx = find(unitsFound == 25002);

showUnitsIdx = idx; %find(unitsFound == maxNeuronID);
for i=1:length(showUnitsIdx)
    uIdx = showUnitsIdx(i);
    uID = unitsFound(uIdx);
    uElGroup = (uID - mod(uID,1000))/1000;
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
save(fullfile(resultDataPath, 'Sorting1.mat'), 'sBM', 'BM', 'SCall', 'EL', 'ELNR');

%% EXPORT ELECTRODES OF 40 SELECTED NEURONS FOR CONFIG GENERATION SCRIPT
elIdx = 8:17;
selNeurons = sBM;
ELs4Michele = [];
ELs4Michele(:,:,1) = reshape(EL(selNeurons(:,elIdx),1), [size(selNeurons,1) length(elIdx)]);
ELs4Michele(:,:,2) = reshape(EL(selNeurons(:,elIdx),2), [size(selNeurons,1) length(elIdx)]);
ELs4Michele(:,:,3) = reshape(ELNR(selNeurons(:,elIdx),1), [size(selNeurons,1) length(elIdx)]);
save(fullfile(resultDataPath, 'Electrodes4Michele.mat'), 'ELs4Michele');

%% ROUTE CONFIG
mysort.mea.routeExperimentConfig(resultDataPath, ELs4Michele);
