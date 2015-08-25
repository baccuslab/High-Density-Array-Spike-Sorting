
function R = compareGDFs(gdf1, gdf2, maxJitter, maxShift)
    % Make sure the GDFs are really sorted according to the timestamp
    gdf1 = sortrows(gdf1,2);
    gdf2 = sortrows(gdf2,2);
    st1 = mysort.spiketrain.fromGdf(gdf1);
    st2 = mysort.spiketrain.fromGdf(gdf2);
    
    R = mysort.spiketrain.compare(st1, st2, maxJitter, maxShift);
end