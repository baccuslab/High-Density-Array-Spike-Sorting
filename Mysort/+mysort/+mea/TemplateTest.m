%% Build fake Data
nChannel = 10;
X = randn(10000, nChannel);
samplesPerSecond = 20000;
name = 'Config1';
epos = randi(1000, [nChannel 2]);
enrs = 1:nChannel;

% Build one Data Source for First Configuration:
DS1 = mysort.ds.Matrix(X, samplesPerSecond, name, epos, enrs);
nSpikes = 1000;
spikeTimes = randi(9000, nSpikes) + 100;
nNeurons = 8;
spikeNeuronIDs = randi(nNeurons, nSpikes);
cutLeft = 10;
cutLength = 40;
gdf = [spikeNeuronIDs(:) spikeTimes(:)];
% Build Spike Sorting Container
SC1 = mysort.spiketrain.SpikeSortingContainer('TestSorting', gdf, ...
        'wfDataSource', DS1,...
        'templateCutLeft', cutLeft,...
        'templateCutLength', cutLength);

% Build a second Data Source for Second Configuration. This configuration
% overlaps with the first one in terms of electrodes
epos(5:end,:) = randi(1000, [6 2]);
enrs(5:end) = enrs(5:end) + nChannel;
DS2 = mysort.ds.Matrix(X, samplesPerSecond, name, epos, enrs);
nSpikes = 1000;
spikeTimes = randi(9000, nSpikes) + 100;
nNeurons = 8;
spikeNeuronIDs = randi(nNeurons, nSpikes);
cutLeft = 10;
cutLength = 40;
gdf = [spikeNeuronIDs(:) spikeTimes(:)];
% Build Second Spike Sorting Container
SC2 = mysort.spiketrain.SpikeSortingContainer('TestSorting', gdf, ...
        'wfDataSource', DS2,...
        'templateCutLeft', cutLeft,...
        'templateCutLength', cutLength);
    
%% Merge Templates
TL1 = SC1.getTemplateList();
TL2 = SC2.getTemplateList();
TLmerged = TL1(1).merge(TL2(1));

mysort.plot.templates2D(TLmerged.waveforms, TLmerged.MultiElectrode.electrodePositions, 20, 0)
hold on
mysort.plot.templates2D(TL1(1).waveforms, TL1(1).MultiElectrode.electrodePositions, 20, 0, 'ah', gca, 'IDs', 2)

TL1(1).MultiElectrode.electrodeNumbers
TL2(1).MultiElectrode.electrodeNumbers
TLmerged.MultiElectrode.electrodeNumbers

