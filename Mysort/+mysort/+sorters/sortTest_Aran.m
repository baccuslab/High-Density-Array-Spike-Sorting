% simple test scenario for the sorter
%% INIT
outPath = '/Users/Aran/Documents/MATLAB/SVN_Matlab_export_R34090/' %directory of your choice
sortingName = 'Test22';
fname = 'testdata_janelia.h5';

%% LOAD
X_unfil = double(hdf5read(fullfile(fname), '/data')');
samplesPerSecond = double(hdf5read(fullfile(fname), '/samples_per_second'));
b  = mysort.mea.filter_design_fir(250, 7000, samplesPerSecond, 101);
X = conv2(X_unfil,b(:),'same');
epos = double(hdf5read(fullfile(fname), '/electrode_pos')');
epos = epos(:,2:3);
enrs = 1:size(epos,1);
gdf = double(hdf5read(fullfile(fname), '/ground_truth_spike_trains')');
gdf(:,2) = samplesPerSecond*gdf(:,2)/10000000;  % convert to samples
name = 'EspenPolytrode4';
DS_unfil = mysort.ds.Matrix(X_unfil, samplesPerSecond, 'unfiltered', epos, enrs);
DS = mysort.ds.Matrix(X, samplesPerSecond, name, epos, enrs);

% Only sort first 8 channels (of 16). Results in speed up for testing and
% 16 is also too much to sort at once
DS.restrictToChannels(1:8);

%% SPIKE SORT
P = struct();
P.spikeDetection.thr = 4.1;
P.spikeCutting.chunkSize = 5000;
P.spikeCutting.blockwise = false;
P.clustering.meanShiftBandWidthFactor = 1.3;
[S P] = mysort.sorters.sort(DS, outPath, sortingName, P);
gdf_sorted = [S.clusteringMerged.ids(:) S.spikeDetectionMerged.ts(:)];

%% EVAL
R = mysort.spiketrain.alignGDFs(gdf, gdf_sorted, 20, 20, 20);
mysort.plot.printEvaluationTable(R);

%% PLOTS
mysort.plot.SliderDataAxes({DS_unfil DS});