%%
idx1 = [1 1 1 1 2 5 1 5 5 3 3 5 4 3]';
idx2 = [2 2 3 2 1 5 2 3 5 4 4 5 1 4]';

[groups occurances idx_occurences] = mysort.util.findCommonClusters([idx1 idx2], 'doPlot', 1)

%%
idx1 = [1 1 2 1 5 3 4 5 1 4 3 1 1]';
idx2 = [2 3 1 2 3 4 1 5 2 1 4 2 2]';
idx3 = [1 2 5 6 1 3 4 1 1 4 4 1 3]';
[groups occurances idx_occurences] = mysort.util.findCommonClusters([idx1 idx2 idx3])


