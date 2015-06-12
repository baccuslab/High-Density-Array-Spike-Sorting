function RES = E13_processing(GT, expsig, save_path) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort the data
    SpDet = mysort.detectors.threshold('method', 'none',...
                                       'threshFactor', 3,...
                                       'cache_path', save_path,...
                                       'cache_name', [GT.source_filename 'Detection'],...
                                       'cacheable', true);
    EM_sorter = mysort.sorters.EMClustering('Tf',           45,...
                                            'cutLeft',      13,...
                                            'nMinCluster',  1,...
                                            'nMaxCluster',  5,...
                                            'EXPSIG',       expsig,...
                                            'nClusterDims', 8,...                                    
                                            'debug',        false,...
                                            'spikeDetector', SpDet);

    gdf = EM_sorter.sort(GT.X);

    RES.EM_sorter = EM_sorter;
    RES.gdf = gdf;
    RES.source_filename = GT.source_filename;