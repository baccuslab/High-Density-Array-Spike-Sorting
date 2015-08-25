function RES = E12_eval(GT, R, maxJitter, maxOverlapDist, maxShift, save_path) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do the evaluation
    gdf1 = GT.gdf;
    gdf2 = R.gdf;
    fprintf('Evaluation of %s..\n', GT.source_filename);
    
    RES = mysort.spiketrain.alignGDFs(gdf1, gdf2, maxJitter, maxShift, maxOverlapDist);
    RES.M = [sum(RES.detErr) sum(RES.detErrO) sum(RES.nCL) ...
             sum(RES.nCLO) sum(RES.totErr) sum(RES.totErrO)];
    RES.source_filename = GT.source_filename;
    fprintf('done.\n');
