%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set path to simulated data files from Rodrigo Quian Quiroga 2004. If you
% downloaded the wave_clus 2.0 package, this path will probably end with:
% /wave-clus-2.0/Simulator/
% The package can be downloaded here:
% http://www.vis.caltech.edu/~rodri/Wave_clus/Wave_clus_home.htm
path_to_simulated_files = E10_quirogaDataPath();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data, use the ground truth information to cut out spikes,
% calculate the noise statistics and the correct templates
% Choose a benchmark and a noise level
benchmarks = {'Easy1', 'Easy2', 'Difficult1', 'Difficult2'};
b = 1; noisestr = '005';
quirogaFile = sprintf('C_%s_noise%s.mat', benchmarks{b}, noisestr);

%print the statistics to that benchmark:
E10_printQuirogaEvaluation(b);

% Set the length of the templates and where they should be cut in respect
% to the spike times stored in the ground truth data
Tf = 71; cutLeft = -10;
GT = E10_preprocessing(path_to_simulated_files, quirogaFile, cutLeft, Tf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort the data using the correct spike templates and noise estimates
GT.C = GT.C + eye(size(GT.C))*GT.C(1,1)/10;
NE = mysort.util.NoiseEstimator(GT.C, Tf);
botm = mysort.sorters.BOTM(NE, Tf, GT.templates, 'upsample', 3, 'spikePrior', .01);

%botm.plotTemplates()
%GT.X = GT.X(:, 3000:3125+200-1);
DH = mysort.ds.Matrix(GT.X', 22000);
gdf = botm.sort(DH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the evaluation
maxJitter = 15;
maxOverlapDist = floor(Tf/2);
maxShift = 0;
R = mysort.spiketrain.alignGDFs(GT.gdf, gdf, maxJitter, maxShift, maxOverlapDist);
fprintf('\nResult for file: %s\n', quirogaFile);
mysort.plot.printEvaluationTable(R);

botm.plotLastChunkSorting('X', GT.X);
% mysort.plot.sortingErrors(R, botm, 'CLO', 'X', GT.X);
figure; mysort.plot.sortingErrors(R, botm, 'CLO');

% botm.plotLastChunkSorting('X', GT.X, 'start', 1422350-50, 'stopp', 1422790-60, 'srate', 1/24000)
% set(gca, 'ylim', [-75 45]);
% xlabel('time [s]');
%%
SSC = mysort.spiketrain.SpikeSortingContainer('Ground truth', GT.gdf, 'templateWfs', GT.XI, 'templateCutLeft', 0);
DH.SpikeSortingContainers = {SSC};
mysort.plot.SliderDataAxes({DH, botm});
