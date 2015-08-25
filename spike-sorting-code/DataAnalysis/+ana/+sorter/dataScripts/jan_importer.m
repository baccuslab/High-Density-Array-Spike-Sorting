%% INIT PATHDEFS
pd = pdefs();

pd.mea1kData
files = {[pd.mea1kData 'jamuelle/131030/data/mostActiveEls/000.raw.h5']
         [pd.mea1kData 'jamuelle/131030/data/mostActiveEls/001.raw.h5']};

nrk_file = {[pd.mea1kData 'jamuelle/131030/config/mostActiveEls/mostActiveEls.mapping.nrk']};

name = 'JAN_131030_mostActiveEls';

inputPath = fullfile(pd.networkTempShare, 'Hillier_2013');

configNTKLists =  {files};
xy_pos         =  {nrk_file};

outfile = fullfile(inputPath, [name '.mat']);

%% SAVE CONFIG FILE FOR SPIKE SORTER
% assert(~exist(outfile, 'file'), 'Does already exist!');
% save(outfile, 'configNTKLists', 'xy_pos');


%% TO START SPIKE SORTING
% go to "startHillerMany2013_grid.m" in the same folder as this script.
% enter the name of the outfile saved in the previous cell in 
% the "expNames" names variable (as a cell with one element)
% The spike sorting deamon must run on one of the submit hosts, to push the
% sorting to the cluster


%% CHECK RESULTS, INIT PATHS
intermPath = fullfile(pd.mea1kIntermediate, 'frankef', name);
resPath = fullfile(intermPath, 'SortOut', 'Sorting1', 'Config1', [name '_proj_1']);
resFile = fullfile(resPath, [name '_proj_1_results.mat']);

res = load(resFile);
h5files = {fullfile(intermPath, 'Preprocessed', '000.raw.h5_prefilter.h5')
           fullfile(intermPath, 'Preprocessed', '001.raw.h5_prefilter.h5')};

%% INIT THE DATA ACESS OBJECT
cmos = mysort.mea.CMOSMEA(h5files);

%% LOAD THE GROUP FILE WHICH SPECIFIES WHICH ELECTRODES WERE SORTED TOGETHER
G = load(fullfile(resPath, 'groupFile.mat'));
nElectrodes = size(cmos, 2);

%% REMOVE NEURONS THAT DO NOT HAVE ENOUGH SPIKES OR THE DESIRED AMPLITUDE
MIN_SPIKES_PER_NEURON = 100;
MIN_AMPLITUDE_PER_NEURON = 30;

gdf = res.gdf_merged;
units = unique(gdf(:,1));
edges = [units(:)-.5 units(:)+.5]';
edges = edges(:);
nSpikesPerNeuron = histc(gdf(:,1), edges);
nSpikesPerNeuron(2:2:end) = [];
min(nSpikesPerNeuron)

ampsPerNeuronAndElectrode = squeeze(max(abs(res.T_merged), [], 1));
ampsPerNeuron = squeeze(max(ampsPerNeuronAndElectrode, [], 1));

removeIdx1 = nSpikesPerNeuron(:) < MIN_SPIKES_PER_NEURON;
fprintf('Removing %d neurons because they had too few spikes\n', length(find(removeIdx1)));


removeIdx2 = ampsPerNeuron(:) < MIN_AMPLITUDE_PER_NEURON;
fprintf('Removing %d neurons because they had too low amplitude\n', length(find(removeIdx2)));

gdf_final = gdf;
removeSpikesIdx = ismember(gdf(:,1), find(removeIdx1 | removeIdx2));
gdf_final(removeSpikesIdx,:) = [];

%% THIS BUILDS THE SPIKE SORTING CONTAINER, WHICH MIGHT BE USEFUL FOR PLOTTING
SSC = mysort.spiketrain.SpikeSortingContainer('Bla', gdf_final, 'wfDataSource', cmos);

%% WITH THESE YOU CAN PLOT SOME SELECTED NEURONS
showUnitIdx = 2500;
showUnitsNames = 344010; %units_final(showUnitIdx);

% AND THIS LOOP ACTUALLY DOES THE PLOTTING
units_final = unique(gdf_final(:,1));
T_final = res.T_merged;
T_final(:,:, removeIdx1 | removeIdx2) = [];
cutLeft = 45;
cutLength = 85; 

for i=1:length(showUnitsNames)
    
    uIdx = find(units_final == showUnitsNames(i));
    uID = units_final(uIdx);
    uElGroup = (uID - mod(uID,1000))/1000;
    
%     nSpikes = find(gdf_final
    spikes = SSC.getWaveforms4UnitIdx(uIdx, cutLeft, cutLength, 200);
    tspikes = mysort.wf.v2t(spikes, nElectrodes);
    EP = cmos.MultiElectrode.electrodePositions;
    groupElIdx = G.groupsidx{uElGroup};
    
    mysort.plot.figure([1200 800])
    ah = axes;
    P.scaling = 0.02;
    P = mysort.plot.waveforms2D(tspikes, EP, 'plotMedian', 1, 'AxesHandle', ah, 'scaling', P.scaling);
    title(sprintf('Unit Name: %d', uID))
    hold on
    plot(EP(groupElIdx,1), EP(groupElIdx,2), 'rx', 'markersize', 14, 'linewidth', 3)    
    plotArgs = {'-', 'color', [.1 .6 .2], 'linewidth', 3};
%     mysort.plot.waveforms2D(T_final(:,:,uIdx), EP, 'plotArgs', plotArgs, 'AxesHandle', ah, 'scaling', P.scaling);   
    
end

%% THIS IS AN ALTERNATIVE WAY TO PLOT THE SPIKES USING THE GUI
WF_all = mysort.wf.BufferedWfManager(cmos, gdf_final(:,2), gdf_final(:,1), [], cutLeft, cutLength);

maxNWaveforms = 50;
WFP_all = mysort.plot.WaveformManagerPlot(WF_all, 'plotMedianColor', [.9 .7 .0],...
                        'plotControls', true, 'maxPlotWaveforms', maxNWaveforms, ...
                        'maxPlotChannels', 6);


