
function R = alignGDFs(gdf1, gdf2, maxJitter, maxShift, maxOverlapDist)
    % For an explanation how to use this function check
    % mysort.spiketrain.align
    %
    % This function is just a wrapper which converts the two gdfs *1) into
    % cell arrays of spike trains and calls mysort.spiketrain.align
    %
    % *1) gdf1 and gdf2, which are matrices of spike trains, first column
    % neuron, second column time of spike in samples and each row
    % represents a spike)
    st1 = {};
    if ~isempty(gdf1)
        gdf1 = double(sortrows(gdf1,2));
        [st1 newids1 ids1] = mysort.spiketrain.fromGdf(gdf1);
    end
    
    st2 = {};
    if ~isempty(gdf2)
        gdf2 = double(sortrows(gdf2,2));
        [st2 newids2 ids2]= mysort.spiketrain.fromGdf(gdf2);
    end
    R = mysort.spiketrain.align(st1, st2, maxJitter, maxShift, maxOverlapDist, ids1, ids2);
end