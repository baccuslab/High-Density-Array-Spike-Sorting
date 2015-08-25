%% INIT 
%matlabpool(12)
IS_OFFLINE = false;

pd = pdefs();
% make sure you have the path of a folder with all .ntk files WITH A SIN% COMMON configuration in the variable "confPath"
% expName = 'Cortex mice';
expName = '150728online'; %'150630'; %'150630offline'; %'150519offline'; %150512wnoi %150512offline % 150519bscreen
% expName = '150630offline'; 
sortingName = 'Sorting1';
expDataPath          = fullfile(pd.mea1kData, 'Antonia', expName);
intermediateDataPath = fullfile(pd.mea1kIntermediate, 'Antonia', expName);
outDataPath          = fullfile(pd.serverData, 'Antonia', expName);

h5Path = [];

toni_create_flist(expName, outDataPath)

INTERMEDIATE_DATA_PATH = fullfile(pd.mea1kIntermediate, 'frankef', 'Antonia');
CONFIG_FILE_LOCATION = outDataPath;
% resultFolder = fullfile(outDataPath, 'Results');
submitForReal = true;
% submitForReal = false;
runLocally = false;
bCorrectForMissingFrames = false;
expNames = {expName};

gridSpikeSorting();


%% EXPORT RESULTS FOR ANTONIA
flist = PATHDEFS(1).h5FileLists{1};
configIdx = 1;

outPath = fullfile(PATHDEFS(1).sortOutPaths{1});
[SSC S compoundMea G] = ana.toni.inspectSorting(PATHDEFS(1).jobName, outPath, flist);   

% Get spike trains and chop into different files/conditions
L = compoundMea.X.getAllSessionsLength();

% break gdf apart
start_end_times = [0 cumsum(L)];
mgdf = mysort.spiketrain.gdf2multiSessionGdf(S.gdf_merged, start_end_times);
R = {};

units = unique(mgdf(:,1));
for k=1:length(L)
    gdf = mgdf(mgdf(:,3)==k,[1:2]);

    R{k} = gdf;
end

for file=1:length(flist)
    [a,b,c] = fileparts(flist{file});
    cutOffIdx = strfind(b, '.');
    assert(~isempty(cutOffIdx), 'problem!');
    cutOffIdx = min(cutOffIdx);
    assert(cutOffIdx>2, 'Problem!')
    fname = [b(1:cutOffIdx-1) '.mat'];
    resultFile = fullfile(outDataPath, fname);
    
    firstFrameNumber = double(hdf5read(flist{file}, '/Sessions/Session0/frame_numbers/first_fn'));

    spikeTrainsInSamples = R{file};
    spikeTrainsInFrameNumbers = spikeTrainsInSamples + firstFrameNumber - 1;
    save(resultFile, 'firstFrameNumber', 'G', 'spikeTrainsInSamples', 'spikeTrainsInFrameNumbers', 'units');
end


if IS_OFFLINE
    disp('OFFLINE analysis done!')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Starting ONLINE analysis!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% EXTRACT FRAMENUMBERS
logFolderNames = {'control', 'PSEM', 'wash'};

clear STIMINFO;
for i=1:length(PATHDEFS(1).h5FileLists)
    for k=1:length(PATHDEFS(1).h5FileLists{i})
        rawh5File = PATHDEFS(1).configNTKLists{i}{k};
        h5File  = PATHDEFS(1).h5FileLists{i}{k};
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
        SI.startSampleTime = [];
        SI.endSampleTime = [];
        for j=1:length(StimInfo.stimulus_frame_info)
            if (strcmp(StimInfo.stimulus_frame_info{j}.stimulus_type, 'show_image') || ...
                strcmp(StimInfo.stimulus_frame_info{j}.stimulus_type, 'white_noise')) && ...
               ~StimInfo.stimulus_frame_info{j}.is_last
               counter   = double(StimInfo.stimulus_frame_info{j}.counter);
               cmd_index = counter*2 + 1;
               framenumber = double(CMD(cmd_index, 1));
               sampleTime  = framenumber - firstFrameNumber + 1;
               
               SI.stimInfoIndex = [SI.stimInfoIndex j];
               SI.counter = [SI.counter counter];
               SI.cmd_index = [SI.cmd_index cmd_index];
               SI.framenumber = [SI.framenumber framenumber];
               SI.startSampleTime = [SI.startSampleTime sampleTime];
               
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

%% Extract Good Electrodes for all Templates
clear GE;
disp('Computing good Electrodes For Each Template And Config...');
tic

configIdx = 1;

unitsFound = SSC.unitNames;
cutLeft = 10;
cutLength = 45; 
nElectrodes = size(compoundMea, 2);

% Get spike trains and chop into different files/conditions
L = compoundMea.X.getAllSessionsLength();
% break gdf apart
start_end_times = [0 cumsum(L)];

R_stimuli = {};
SC = [];
units = unique(mgdf(:,1));
for k=1:length(L)
    gdf = mgdf(mgdf(:,3)==k,[1:2]);
    for u=1:length(units)
        SC(k,u) = length(find(gdf(:,1)==units(u)));
    end
    TS = [STIMINFO(k,configIdx).startSampleTime(:) STIMINFO(k,configIdx).endSampleTime(:)];
    R_stimuli{k} = mysort.spiketrain.gdf2multiSessionGdf(gdf, TS, [], 0);
end

%% Get Electrode Numbers and best Electrodes for Templates 
elNr = compoundMea.MultiElectrode.electrodeNumbers;
T = S.T_merged;
TemplateAmplitudes = squeeze(max(max(abs(T), [], 1), [], 2));

[mi ma] = mysort.wf.tMinMaxPerTemplate(T);
mi = abs(mi); ma = abs(ma);
M = max(mi, ma);
SortedTemplateAmplitudeIndex = zeros(size(M));
SortedTemplateAmplitudes = SortedTemplateAmplitudeIndex;
for i=1:size(SortedTemplateAmplitudeIndex,1)
    [SortedTemplateAmplitudes(i,:) SortedTemplateAmplitudeIndex(i,:)] = sort(M(i,:), 'descend');
end
%% COMPUTE VARIABILITY OF RESPONSES OVER TRIALS
CONTROL = 1;
PSEM = 2;
WASH = 3;
FIDX_SET = [CONTROL PSEM WASH];
nStimulusRepetitions = length(STIMINFO(CONTROL).stimInfoIndex);
StimulusTime  = (STIMINFO(CONTROL).endSampleTime(1)-STIMINFO(CONTROL).startSampleTime(1)) / 20000;
nConditions = length(FIDX_SET);


window = .15;
sampling = .05;
timePeriod = [0 StimulusTime]; % Trial length in seconds
nTimeSamples = ceil(((timePeriod(2)-timePeriod(1)) / sampling))+1;
FR_per_trial_unit_condition = zeros(nTimeSamples, nStimulusRepetitions, length(units), nConditions);

for file = 1:nConditions
    for u = 1:length(units)
        uID = units(u);
        for k = 1:nStimulusRepetitions
            stidx = R_stimuli{file}(:,3) == k & R_stimuli{file}(:,1) == uID;
            st = R_stimuli{file}(stidx,2);
            [fr_, x] = mysort.spiketrain.toFiringRate(st/20000, window, sampling, timePeriod);
            FR_per_trial_unit_condition(:, k, u, file) = fr_;
        end
    end
end

%% COMPUTE MODULATION OF SPIKE TRAINS
FF_per_time_neuron_condition  = zeros(nTimeSamples, length(units), nConditions);
FF_of_mean_neuron_condition   = zeros(length(units), nConditions);
for file = 1:nConditions
    for u = 1:length(units)
        % COMPUTE FANO FACOR FOR EACH CONDITION
        meanFR = squeeze(mean(FR_per_trial_unit_condition(:, :, u, file), 2));
        varFR  = squeeze(var(FR_per_trial_unit_condition(:, :, u, file), [], 2));
        ff    = varFR./meanFR;
        FF_per_time_neuron_condition(:, u, file) = ff;
        FF_of_mean_neuron_condition(u, file) = var(meanFR)/mean(meanFR);
    end
end

if 0
    %% FIRING RATES
    cols = {'b', 'r', 'g'};
    for uidx = 10:20;
        figure
        subplot(2,1,1)
        hold on
        for i=1:3
            plot(x, squeeze(FR_per_trial_unit_condition(:, :, uidx, i)), 'color', cols{i});
        end
%         title(['FF of mean: ' sprintf(' %.2f', FF_of_mean_neuron_condition(uidx,:))])

        subplot(2,1,2)
        hold on
        for i=1:3
%             plot(x, FF_per_time_neuron_condition(:, uidx, i), 'color', cols{i});
        end
    end 
    %% FANO FACTORS
    cols = {'b', 'r', 'g'};
    for uidx = 10:20;
        figure
        subplot(2,1,1)
        hold on
        for i=1:3
            plot(x, squeeze(FR_per_trial_unit_condition(:, :, uidx, i)), 'color', cols{i});
        end
        title(['FF of mean: ' sprintf(' %.2f', FF_of_mean_neuron_condition(uidx,:))])

        subplot(2,1,2)
        hold on
        for i=1:3
            plot(x, FF_per_time_neuron_condition(:, uidx, i), 'color', cols{i});
        end
    end
end


%% TODO: ONLY DO THIS FOR PSEM WASH & CONTROL
SC_ = SC+1;
Activity         = sum(SC([CONTROL PSEM WASH],:))/nStimulusRepetitions/3;
SilenceIndex   = (SC_(PSEM,:)+1)./(mean(SC_([CONTROL WASH],:))+1);
StabilityIndex = (SC_(CONTROL,:)+1)./(SC_(WASH,:)+1);

ResponseIndex  = min(FF_of_mean_neuron_condition(:,[CONTROL WASH]), [], 2);

g = struct();
g.elNr = elNr;
g.unitsFound = unitsFound;
g.elPos = compoundMea.MultiElectrode.electrodePositions;
g.elLabels = compoundMea.MultiElectrode.electrodeLabels;

g.goodElectrodes = SortedTemplateAmplitudeIndex;
g.goodElectrodesAmplitudes = SortedTemplateAmplitudes;
g.goodElectrodesAmplitudesMax = max(abs(g.goodElectrodesAmplitudes), [], 2);
g.mgdf = mgdf;
g.SingleFileGDF = R;
g.FileLength = L;
g.Units = units;
g.SpikeCountsPerFile = SC;
g.SilenceIndex = SilenceIndex;
g.StabilityIndex = StabilityIndex;
g.Activity  = Activity;
g.TemplateAmplitudes = TemplateAmplitudes;
g.SingleStimulusGDF = R_stimuli;
g.ResponseIndex = ResponseIndex;

GE = g;
clear g


if 0
    %%
    mysort.plot.figure([1400 700]);
    subplot(2,3,1)
    hist(GE.Activity, 100);
    title('Activity');
    
    subplot(2,3,2)
    hist(GE.SilenceIndex, 100);
    title('SilenceIndex');
    
    subplot(2,3,3)
    hist(GE.StabilityIndex, 100);
    title('StabilityIndex');
    
    subplot(2,3,4)
    hist(GE.TemplateAmplitudes, 100);
    title('TemplateAmplitudes');
    
    subplot(2,3,5)
    hist(GE.ResponseIndex, 100);
    title('FF of mean');
end   


%% PRESELECTION OF NEURONS BASED ON TEMPLATE AND ACTIVITY
thresholds = {'Activity',            '>', 10
              'TemplateAmplitudes',  '>', 25
              'SilenceIndex',        '<', 100
              'ResponseIndex',       '>', 1
              'StabilityIndex',      '>', .5
              'StabilityIndex',      '<', 1.5};

nThr = size(thresholds,1);
indices = zeros(length(GE.Units), nThr);
for ti = 1:nThr
    if strcmp(thresholds{ti, 2}, '>')
        indices(:, ti) = GE.(thresholds{ti,1}) > thresholds{ti,3};
    else
        indices(:, ti) = GE.(thresholds{ti,1}) < thresholds{ti,3};
    end
end

%
PRESELECTION_INDEX = find(all(indices, 2));

or_idx  = sum(indices==0,1);
xor_idx = sum(indices==0 & repmat(sum(indices==0,2) == 1, [1, nThr]));

mysort.plot.printTable([or_idx; xor_idx], 'colLabel', thresholds(:,1), 'rowLabel', {'or', 'xor'});
fprintf('##### %d Neurons preselected\n', length(PRESELECTION_INDEX));

%% COMPUTE RECEPTIVE FIELDS FROM WHITE NOISE
VISUAL_FRAMERATE = 30; % Hz
WHITENOISE_INDEX = 4;
HIST_LENGTH = 11;
rawh5File = PATHDEFS.configNTKLists{configIdx}{WHITENOISE_INDEX};
stimFile = [rawh5File(1:end-6) 'mat'];
StimInfo = load(stimFile);    

gdf = R_stimuli{WHITENOISE_INDEX};

% throw away spikes that did not happen during the stimulus
gdf(gdf(:,3)==0,:) = [];
stimStartTime = STIMINFO(WHITENOISE_INDEX, configIdx).startSampleTime;
stimEndTime = STIMINFO(WHITENOISE_INDEX, configIdx).endSampleTime;
stimLenInSamples = stimEndTime - stimStartTime;
stimLenInSec = stimLenInSamples/20000;

% WNST = StimInfo.stimulus_frame_info{1, STIMINFO(WHITENOISE_INDEX, configIdx).stimInfoIndex}.parameters.color;
duration_in_theory = StimInfo.stimulus_frame_info{1, 1}.parameters.duration;
dd = load(fullfile(pd.serverData, 'Antonia', 'white noise stim info', '75um', 'fragment_WhiteNoiseParameters_00001_0_colors.mat'));
nFrames = round(duration_in_theory*60);
WNST = dd.colors(1:nFrames,:,:);

% THIS IS A HACK !!   SHOULD BE GIVEN BY SETUP, BUT IS UNCLEAR AT THE MOMENT
VISUAL_FRAMERATE = nFrames/stimLenInSec;

stimLenInFrames = stimLenInSec*VISUAL_FRAMERATE;
assert(nFrames >= stimLenInFrames, 'Wrong number of Stimulus Frames');


U = GE.Units;

nPreselectedNeurons = length(PRESELECTION_INDEX);
RF = zeros(size(WNST,2), size(WNST,3), HIST_LENGTH, nPreselectedNeurons);
sRF = RF;
K = [.1  .5  .1
     .5  1.5  .5
     .1  .5  .1];
K = K/norm(K(:));

for uips = 1:nPreselectedNeurons
    ui = PRESELECTION_INDEX(uips);
    st = gdf(gdf(:,1) == U(ui),2);
    if isempty(st)
        fprintf('NO spike found for unit id:%d idx:%d\n', U(ui), ui);
        continue
    end
    % compute spike times as index into movie frame
    st = round(st/20000 * VISUAL_FRAMERATE);

    for f = 1:HIST_LENGTH
        stidx = st-(f-1);
        stidx(stidx<1) = [];
        RF(:,:,f,uips) = mean(WNST(stidx,:,:), 1)-.5;
        sRF(:,:,f,uips) = conv2(RF(:,:,f,uips), K, 'same');
    end
end  
%%
RFcenters = [];
CellTypes = {};
for uips = 1:nPreselectedNeurons
    RFmax = max(sRF(:,:,:, uips), [], 3);
    RFmin = min(sRF(:,:,:, uips), [], 3);
    [a_max,b_max, val_max] = mysort.util.matrixArgMax(RFmax);
    [a_min,b_min, val_min] = mysort.util.matrixArgMax(-RFmin); 
    val_min = -val_min;
    t_max = find(sRF(a_max,b_max,:,uips) == val_max);
    t_min = find(sRF(a_min,b_min,:,uips) == val_min);
    
    [maxAmp, maxFrame] = max(squeeze(abs(sRF(a,b,:,uips))));
    val = sRF(a,b, maxFrame, uips);
    [centroid, type] = ana.antonia.classify2Dmovie(squeeze(RF(:,:,:,uips)), [], maxFrame);
    
    myRf = sRF(:,:,:,uips);
    bHasOn  = val_max >  6.5*std(myRf(:));
    bHasOff = val_min < -6.5*std(myRf(:));
    
    celltype = '';
    if ~bHasOn && ~bHasOff && abs(val_max) > abs(val_min)
        celltype = 'ON';
    elseif ~bHasOn && ~bHasOff && abs(val_max) > abs(val_min)
        celltype = 'OFF';
    elseif bHasOn && ~bHasOff
        celltype = 'ON';
    elseif ~bHasOn && bHasOff
        celltype = 'OFF';
    elseif t_max > t_min
        celltype = 'ON';
    else
        celltype = 'OFF';
    end
    
    if strcmp(celltype, 'ON')
        rfcenter = [a_max b_max  1 ];
    elseif strcmp(celltype, 'OFF')
        rfcenter = [a_min b_min  0 ];
    else
        rfcenter = [a_min b_min  -1 ];
    end
    
    CellTypes{uips} = celltype;    
    
    RFcenters(uips,:) = rfcenter;
end
save(fullfile(outDataPath, 'rf_centers.mat'), 'RFcenters', 'PRESELECTION_INDEX', 'CellTypes');
disp('Done ONLINE!');
        
%%
if 0 
    %% PLOT RFs !
    for uidx = size(sRF,4)-1
        RFmax = max(sRF(:,:,:, uidx), [], 3);
        RFmin = min(sRF(:,:,:, uidx), [], 3);
        P = ana.antonia.plotRF(sRF(:,:,:,uidx));
        imagesc(RFmin, 'parent', P.ah(end-1));
        hold on
        plot(P.ah(end-1), RFcenters(uidx,2), RFcenters(uidx,1), 'k.', 'markersize', 22)
        plot(P.ah(end-1), RFcenters(uidx,2), RFcenters(uidx,1), 'w.', 'markersize', 19)
        imagesc(RFmax, 'parent', P.ah(end));
        hold on
        plot(P.ah(end), RFcenters(uidx,2), RFcenters(uidx,1), 'k.', 'markersize', 22)
        plot(P.ah(end), RFcenters(uidx,2), RFcenters(uidx,1), 'w.', 'markersize', 19)
                
        mi = min(RFmin(:));
        ma = max(RFmax(:));
        set(P.ah, 'clim', [mi ma]);
        
        mysort.plot.figureTitle(CellTypes{uidx});
    end
end

if 0
    %% PLOT FIRING RATES
    cols = {'b', 'r', 'g'};
    for ui = size(sRF,4)-1
        uidx = PRESELECTION_INDEX(ui);
        figure
        subplot(2,1,1)
        hold on
        for i=1:3
            plot(x, squeeze(FR_per_trial_unit_condition(:, :, uidx, i)), 'color', cols{i});
        end
%         title(['FF of mean: ' sprintf(' %.2f', FF_of_mean_neuron_condition(uidx,:))])

        subplot(2,1,2)
        hold on
        for i=1:3
%             plot(x, FF_per_time_neuron_condition(:, uidx, i), 'color', cols{i});
        end
    end     
end
return











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
    str = sprintf('exp: %s, cellidx: %d, cell unit: %d', expName, sFINAL_SELECTION_INDEX(i), sBM(i,2));
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


