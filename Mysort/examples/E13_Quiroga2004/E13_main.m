%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING
% This script will run for a while depending on your computer speed.
% Furthermore, it will produce a lot of data on your hard disk
error('This file is out of date and needs to be fixed. The dataHandle class should not be used anymore, instead the ds.DataSourceInterface!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INIT
addpath(fullfile('..','E10_Quiroga2004'), fullfile('..','E11_Quiroga2004'),...
        fullfile('..','E12_Quiroga2004'));
path_to_simulated_files = E10_quirogaDataPath();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING
quirogaFiles = E11_getAllQuirogaFiles(path_to_simulated_files);

% Set the length of the templates and where they should be cut in respect
% to the spike times stored in the ground truth data
Tf = 71; cutLeft = -10;
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

exp_sig = [1 2 3 4 5];
parameter_names = {'Expected Std'};
parameter_lists = {exp_sig};

processing = mysort.util.AnalysisHandler('ProcessingE13', @E13_processing, ...
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

evalhandler = mysort.util.AnalysisHandler('EvalE13', @E12_eval, ...
                parameter_names, parameter_lists, processing);
evalhandler.start();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT
evalhandler.surf('M', {'Expected Std', 'Filenames', 5});
zlabel('Total Errors');

evalhandler.surf('totErrO', {'Expected Std', 'Filenames', 1},...
                 'wrapFunc', @sum);

evalhandler.surf('nFP', {'Expected Std', 'Filenames', 1},...
                 'wrapFunc', @sum);
             
evalhandler.surf('detErr', {'Expected Std', 'Filenames', 1},...
                 'wrapFunc', @sum);
             
evalhandler.surf('nSt2', {'Expected Std', 'Filenames', 1});

evalhandler.surf('detErr', {'Expected Std', 'Filenames', 1});

evalhandler.surf('nCL', {'Expected Std', 'Filenames', 1});

evalhandler.surf('nTP', {'Expected Std', 'Filenames', 1},...
                 'wrapFunc', @sum);


