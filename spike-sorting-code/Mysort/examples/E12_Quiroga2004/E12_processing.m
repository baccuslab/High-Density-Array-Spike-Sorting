function RES = E12_processing(GT, spikePrior, condnr, diagonalLoading, save_path) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort the data
    Tf = size(GT.templates, 2);
    NE = mysort.util.NoiseEstimator(GT.C, Tf);
    botm = mysort.sorters.BOTM(NE, Tf, GT.templates, 'upsample', 3,...
                               'spikePrior', spikePrior, ...
                               'minCondNumber', condnr,...
                               'diagonalLoading', diagonalLoading);
    
    gdf = botm.sort(GT.X);
    RES.botm = botm;
    RES.gdf = gdf;
    RES.source_filename = GT.source_filename;