addpath(fullfile('..','E10_Quiroga2004'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set path to simulated data files from Rodrigo Quian Quiroga 2004. If you
% downloaded the wave_clus 2.0 package, this path will probably end with:
% /wave-clus-2.0/Simulator/
path_to_simulated_files = E10_quirogaDataPath();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data, use the ground truth information to cut out spikes,
% calculate the noise statistics and the correct templates
quirogaFiles = E11_getAllQuirogaFiles(path_to_simulated_files);

% Set the length of the templates and where they should be cut in respect
% to the spike times stored in the ground truth data
Tf = 71; cutLeft = -10;

upsample = 3;

maxJitter = 15;
maxOverlapDist = floor(Tf/2);
maxShift = 0;
M = []; R ={};
for i=1:length(quirogaFiles)
    quirogaFile = quirogaFiles{i};
    GT = E10_preprocessing(path_to_simulated_files, quirogaFile, cutLeft, Tf);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort the data

    Covest.CCol = mysort.noise.Cte2Ccol(GT.C, 1);
    alpha = .5;
    Covest.CCol = mysort.noise.Cte2Ccol((1-alpha)*GT.C + alpha*diag(diag(GT.C)), 1);
%     NE = mysort.util.NoiseEstimator(GT.C, Tf);
%     botm = mysort.sorters.BOTM(Covest, Tf, GT.templates, 'upsample', 3, 'spikePrior', .01, 'chunk_size', 500000);
    botm = mysort.sorters.BOTM(Covest, Tf, GT.templates, 'upsample', 1, 'spikePrior', .0001, 'useSIC', false, 'minEpochLength', 8, 'chunk_size', 500000);
    gdf = botm.sort(GT.X);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do the evaluation
    R{i} = mysort.spiketrain.alignGDFs(GT.gdf, gdf, maxJitter, maxShift, maxOverlapDist);
    M(i,:) = [sum(R{i}.detErr) sum(R{i}.detErrO) sum(R{i}.nCL) sum(R{i}.nCLO) sum(R{i}.totErr) sum(R{i}.totErrO)];
    if 0
        DSGT = mysort.ds.Matrix(GT.X', 24000, 'Ground Truth');
        SSCGT = mysort.spiketrain.SpikeSortingContainer('GT', GT.gdf, 'wfDataSource', DSGT);
        DSGT.addSpikeSorting(SSCGT);
        mysort.plot.SliderDataAxes({DSGT, botm});
        mysort.plot.sortingErrors(R{i}, botm, 'CL');
    end    
end

[Q bigTable] = E10_QuirogaPerformanceIn2004Paper();

% Resort the rows in M to match those of the big table
M = M([9:end 1:8],:); rowLabel = quirogaFiles([9:end 1:8]);
bigTable = [bigTable M];

% Remove the Easy1 noise 025 to 04. We dont have the detection errors for 
% the wave_clus analysis and cant compare.
smallTable = bigTable([1:4 9:end],:); rowLabel2 = rowLabel([1:4 9:end]);

% Print a nice table that contrasts the Detection Errors (Det) the 
% Classification Errors (Cla) and the Total Error number (Tot) for 
% the Wave_clus analysis (Q) and for the botm (B)
colLabel = {'QDet', 'QCla', 'QTot', 'BDet', 'BCla', 'BTot'};                         
sumTableSmall = [sum(smallTable(:,3:5),2) smallTable(:,7) sum(smallTable(:,[3:5 7]),2) smallTable(:,[8 10 12])];
sumTableBig = [sum(bigTable(:,3:5),2) bigTable(:,7) sum(bigTable(:,[3:5 7]),2) bigTable(:,[8 10 12])];
mysort.plot.printTable(sumTableSmall, 'rowLabel', rowLabel2, 'colLabel', colLabel,...
                                      'printColSum', 1, 'hlineAfterRows', [4 8 12],...
                                      'vlineAfterCols', 3);                
mysort.util.csvwrite(sprintf('results_small_up%d.csv', upsample), sumTableSmall);

mysort.plot.printTable(sumTableBig, 'rowLabel', rowLabel, 'colLabel', colLabel,...
                                    'printColSum', 1, 'hlineAfterRows', [8 12 16],...
                                    'vlineAfterCols', 3);
mysort.util.csvwrite(sprintf('results_up%d.csv', upsample), sumTableBig);
