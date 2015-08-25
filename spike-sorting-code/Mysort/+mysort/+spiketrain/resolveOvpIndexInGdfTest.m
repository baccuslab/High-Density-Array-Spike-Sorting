gdf = [1 100];
ovpIdx = [1 2 3 -5];
    %   this means that spikes of unit 1 should be replaced by two spikes
    %   each, one of unit 2 with the same time point as the original spike
    %   of unit 1, and one spike of unit 3 with a shift of 5 samples
    %   leftwards
[gdfresolved bResolved] = mysort.spiketrain.resolveOvpIndexInGdf(gdf, ovpIdx)

    %   gdfresolved = [3 95
    %                  2 100];
    
    
gdf = [1 100
       2 200
       3 300
       4 400
       5 500
       3 600
       4 700
       5 800
       6 900
      12 1000
       7 1100];
ovpIdx = [5 1 2 -1
          6 1 2  0
          7 1 2  1
          8 1 3 -1
          9 1 3  0
         10 1 3  1
         11 2 3 -1
         12 2 3  0
         13 2 3  1];
[gdfresolved bResolved] = mysort.spiketrain.resolveOvpIndexInGdf(gdf, ovpIdx)   
