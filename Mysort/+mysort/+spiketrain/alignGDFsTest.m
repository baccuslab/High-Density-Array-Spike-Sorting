% function alignGDFsTest()
    gdf1 = [1 100
            1 120
            1 140
            1 160
            1 200
            1 500
            1 1100
            
            
            2 60
            2 300
            2 600
            
            3 400
            3 800];

    gdf2 = [1 100
            1 120
            1 140
            2 160
            1 200
            1 500
            1 1100
            
            2 60
            2 300
            2 600
            2 1100
            
            3 400
            3 800];

    R = mysort.spiketrain.alignGDFs(gdf1, gdf2, 10, 10, 10);
    mysort.plot.printEvaluationTable(R, 'ignoreEmptySorted', 0);
    mysort.plot.alignment(R)
    
    R.similaritySt1St2_normalized
    R.similarityBeforeAssignmentSt1St2_normalized
    
    R.similaritySt1St2
    R.similarityBeforeAssignmentSt1St2
    disp('Note the difference between R.similaritySt1St2 and R.similarityBeforeAssignmentSt1St2 caused by the spike at 1100');
    disp('This difference comes from the fact that during the assignment this spike is removed from the comparison of St1_1 to St2_2,');
    disp('since this spike is already assigned to a spike in St2_1');
    
    figure;
    imagesc(R.similarityBeforeAssignmentSt1St2_normalized);
    hc = colorbar; ylabel(hc, 'Matching Ratio');
    xlabel('Sorted Spike Trains');
    ylabel('Ground Truth Spike Trains');