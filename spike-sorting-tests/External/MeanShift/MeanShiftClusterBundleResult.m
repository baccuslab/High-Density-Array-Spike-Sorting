function [ids clustCent] = MeanShiftClusterBundleResult(X, cluster2dataCell, minSpikes)
    nPerClus = cellfun(@length, cluster2dataCell);
    fprintf('Clusters with at least %d spikes: %d ;  ', minSpikes, sum(nPerClus>=minSpikes));
    numClust = length(cluster2dataCell);
    ids = zeros(size(X,1),1);
    nClass = 1;
    for k = 1:min(numClust)
        myMembers = cluster2dataCell{k};
        if length(myMembers) > minSpikes
            ids(myMembers) = nClass;
            nClass = nClass+1;            
        end
    end
    fprintf('Spikes in Cluster 0: %d \n', sum(ids==0));
    clustCent = mysort.util.calculateClassMeans(X, ids);