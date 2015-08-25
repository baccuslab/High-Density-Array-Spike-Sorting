
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING
% This script will run for a while depending on your computer.
% Furthermore, it will produce a lot of data on your hard disk
error('This file is out of date and needs some serious rework. The way the mysort.sorters.BOTM class is used is different now!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INIT
addpath(fullfile('..','E10_Quiroga2004'), fullfile('..','E11_Quiroga2004'));
path_to_simulated_files = E10_quirogaDataPath();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING
quirogaFiles = E11_getAllQuirogaFiles(path_to_simulated_files);

% Set the length of the templates and where they should be cut in respect
% to the spike times stored in the ground truth data
Tf = 71; cutLeft = -8;
parameter_names = {'FilePath', 'Filenames', 'cutLeft', 'Tf'};
parameter_lists = {{path_to_simulated_files}, quirogaFiles, cutLeft, Tf};
savePath = 'C:\Data\AnalyseQuirogaTotal\';
preprocessing = mysort.util.AnalysisHandler('Preprocessing', @E10_preprocessing, ...
                            parameter_names, parameter_lists, savePath, ...
                            'resultsInSingleFiles', 1,...
                            'keepResultsInMemory', 1, ...
                            'ignoreParameterForLoadControl', [1]);
preprocessing.start();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESSING
conds = [10000 5000 1000 500 200 100 50 10 5];
spikePrior = [.00001 .0001 .001 0.01 .1];
diagonalLoading = {'DL'};
parameter_names = {'Spike Prior', 'min. Condition Number', 'diagonalLoading'};
parameter_lists = {spikePrior, conds, diagonalLoading};

processing = mysort.util.AnalysisHandler('ProcessingE12', @E12_processing, ...
    parameter_names, parameter_lists, preprocessing, ...
                            'resultsInSingleFiles', 1,...
                            'keepResultsInMemory', 1);
processing.start();



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATION
maxJitter = 15;
maxOverlapDist = floor(Tf/2);
maxShift = 0;
parameter_names = {'max Jitter', 'max Overlap Distance', 'max Shift'};
parameter_lists = {maxJitter, maxOverlapDist, maxShift};

evalhandler = mysort.util.AnalysisHandler('EvalE12', @E12_eval, ...
                parameter_names, parameter_lists, processing);
evalhandler.start();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT
evalhandler.surf('M', { 'min. Condition Number', 'Filenames', 5},...
                 'constraints', {'diagonalLoading', 'DL', 'Spike Prior', .01});
zlabel('Total Errors');