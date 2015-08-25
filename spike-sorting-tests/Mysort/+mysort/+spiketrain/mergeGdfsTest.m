gdf1 = [1 100
        2 110
        1 150
        4 500];
gdf2 = [3 200
        1 400
        5 600];
    
mgdf = mysort.spiketrain.mergeGdfs({gdf1, gdf2})