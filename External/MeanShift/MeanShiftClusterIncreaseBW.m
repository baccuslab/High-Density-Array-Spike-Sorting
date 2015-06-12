function [clustCent,data2cluster,cluster2dataCell] = MeanShiftClusterIncreaseBW(dataPts,bandWidth,plotFlag, nMinPointsInCluster, maxNBandwidthIncreases, bandwidthIncreaseFactor, maxNGradientSpikes)
%perform MeanShift Clustering of data using a flat kernel
%
% ---INPUT---
% dataPts           - input data, (numDim x numPts)
% bandWidth         - is bandwidth parameter (scalar)
% plotFlag          - display output if 2 or 3 D    (logical)
% ---OUTPUT---
% clustCent         - is locations of cluster centers (numDim x numClust)
% data2cluster      - for every data point which cluster it belongs to (numPts)
% cluster2dataCell  - for every cluster which points are in it (numClust)
% 
% Bryan Feldman 02/24/06
% MeanShift first appears in
% K. Funkunaga and L.D. Hosteler, "The Estimation of the Gradient of a
% Density Function, with Applications in Pattern Recognition"
%
% EDIT: Felix Franke 01.07.2013
% If points end up in their own cluster, increase bandwith iteratively to
% also assign those points to bigger clusters
% Parameter:
%   nMinPointsInCluster - this is the criterium when to stop increasing the
%                         bandwidth. If all points are in clusters of at
%                         least this size


    %*** Check input ****
    if nargin < 2
        error('no bandwidth specified')
    end

    if nargin < 3
        plotFlag = false;
    end
    
    if nargin < 7
        maxNGradientSpikes = size(dataPts,2);
    end

    %**** Initialize stuff ***
    [numDim,numPts] = size(dataPts);
    numClust        = 0;
    bandSq          = bandWidth^2;
    initPtInds      = 1:numPts;
    maxPos          = max(dataPts,[],2);                          %biggest size in each dimension
    minPos          = min(dataPts,[],2);                          %smallest size in each dimension
    boundBox        = maxPos-minPos;                        %bounding box size
    sizeSpace       = norm(boundBox);                       %indicator of size of data space
    stopThresh      = 1e-3*bandWidth;                       %when mean has converged
    clustCent       = [];                                   %center of clust
    beenVisitedFlag = zeros(1,numPts,'uint8');              %track if a points been seen already
    numInitPts      = numPts;                               %number of points to posibaly use as initilization points
    clusterVotes    = zeros(1,numPts,'uint16');             %used to resolve conflicts on cluster membership
    clusterWeights  = [];
    
    if maxNGradientSpikes<numPts
        % restrict the estimation of the gradient to a subset of all
        % points. This reduces accurracy but increases speed!
        
    end
    
    bIncreaseBandwidth = true;
    bwIncreaseCounter = 1;
    while bIncreaseBandwidth

        msIteration();
        [val, data2cluster] = max(clusterVotes,[],1);                %a point belongs to the cluster with the most votes
        
        % Check if we have points in clusters that are not big enough. If
        % yes, set those points as unvisited, increase BW and continue
        % clustering but only for those points.
        removeCluster = zeros(1, numClust);
        cluster2dataCell = cell(numClust,1);
        for cN = 1:numClust
            myMembers = find(data2cluster == cN);
            cluster2dataCell{cN} = myMembers;
            if length(myMembers) < nMinPointsInCluster
                removeCluster(cN) = 1;
            end
        end
        if bwIncreaseCounter > maxNBandwidthIncreases || ~any(removeCluster)
            bIncreaseBandwidth = false;
        else
            removeIdx = find(removeCluster);
            clustCent(:,removeIdx) = [];
            clusterVotes(removeIdx,:) = [];
            clusterWeights(removeIdx) = [];
            numClust = numClust - length(removeIdx);
            for i=1:length(removeIdx)
                rIdx = removeIdx(i);
                rIdxPoints = cluster2dataCell{rIdx};               
                beenVisitedFlag(rIdxPoints) = 0;
            end
            initPtInds = find(beenVisitedFlag == 0);
            numInitPts = length(initPtInds);
            bandWidth = bandWidth*bandwidthIncreaseFactor;
            bandSq    = bandWidth^2;
            bwIncreaseCounter = bwIncreaseCounter+1;
        end
    end
    
    

%     %*** If they want the cluster2data cell find it for them
%     if nargout > 2
%         cluster2dataCell = cell(numClust,1);
%         for cN = 1:numClust
%             myMembers = find(data2cluster == cN);
%             cluster2dataCell{cN} = myMembers;
%         end
%     end

    %----------------------------------------------------------------------
    function msIteration()
        numInitPts_old = numInitPts;
        successive_small_clusters = 0;
        t = tic;
        while numInitPts
            tt = toc(t);
            fprintf('Number of active points in set: %d Time since last iteration: %.3f sec\n', numInitPts, tt);
            t = tic();
            tempInd         = ceil( (numInitPts-1e-6)*rand);        % pick a random seed point
            stInd           = initPtInds(tempInd);                  % use this point as start of mean
            myMean          = dataPts(:,stInd);                     % intilize mean to this points location
            myMembers       = [];                                   % points that will get added to this cluster                          
            thisClusterVotes    = zeros(1,numPts,'uint16');         % used to resolve conflicts on cluster membership

            while 1     %loop untill convergence

                sqDistToAll = sum((repmat(myMean,1,numPts) - dataPts).^2);    % dist squared from mean to all points still active
                inInds      = find(sqDistToAll < bandSq);               % points within bandWidth
                thisClusterVotes(inInds) = thisClusterVotes(inInds)+1;  % add a vote for all the points belonging to this cluster


                myOldMean   = myMean;                                   % save the old mean
                myMean      = mean(dataPts(:,inInds),2);                % compute the new mean
                myWeight    = length(inInds);                           % store number of potentially assigned points
                myMembers   = [myMembers inInds];                       % add any point within bandWidth to the cluster
                beenVisitedFlag(myMembers) = 1;                         % mark that these points have been visited

                %*** plot stuff ****
                if plotFlag
                    figure(12345),clf,hold on
                    if numDim == 2
                        plot(dataPts(1,:),dataPts(2,:),'.')
                        plot(dataPts(1,myMembers),dataPts(2,myMembers),'ys')
                        plot(myMean(1),myMean(2),'go')
                        plot(myOldMean(1),myOldMean(2),'rd')
                        pause
                    end
                end

                %**** if mean doesn't move much stop this cluster ***
                if norm(myMean-myOldMean) < stopThresh

                    %check for merge posibilities
                    mergeWith = 0;
                    for cN = 1:numClust
                        distToOther = norm(myMean-clustCent(:,cN));     %distance from posible new clust max to old clust max
                        if distToOther < bandWidth/2                    %if its within bandwidth/2 merge new and old
                            mergeWith = cN;
                            break;
                        end
                    end


                    if mergeWith > 0    % something to merge
                        % find relative weight of cluster means
                        totalWeight = clusterWeights(mergeWith) + myWeight;
                        alpha = clusterWeights(mergeWith)/totalWeight;
                        beta  = myWeight/totalWeight;
                        clustCent(:,mergeWith)       = beta*myMean + alpha*clustCent(:,mergeWith);             %record the max as the mean of the two merged (I know biased twoards new ones)
                        %clustMembsCell{mergeWith}    = unique([clustMembsCell{mergeWith} myMembers]);   %record which points inside 
                        clusterVotes(mergeWith,:)    = clusterVotes(mergeWith,:) + thisClusterVotes;    %add these votes to the merged cluster
                    else    %its a new cluster
                        numClust                    = numClust+1;                   %increment clusters
                        clustCent(:,numClust)       = myMean;                       %record the mean  
                        clusterWeights(1, numClust) = myWeight;
                        %clustMembsCell{numClust}    = myMembers;                    %store my members
                        clusterVotes(numClust,:)    = thisClusterVotes;
                    end

                    break;
                end
            end    
            initPtInds      = find(beenVisitedFlag == 0);           % we can initialize with any of the points not yet visited
            numInitPts      = length(initPtInds);                   % number of active points in set
            if successive_small_clusters > -1
                % Only do this at the very beginning
                if numInitPts_old - numInitPts == 1
                    successive_small_clusters = successive_small_clusters+1;
                else
                    % If one fails, we never look again
                    successive_small_clusters = -1;                
                end
                if successive_small_clusters >= 10 && length(beenVisitedFlag) > 2000;
                    % This is a safety mechanism. If the bandwidth is
                    % completely out of range (too small) each point will fall
                    % in its own cluster. This might take ages though, so stop
                    % here
    %                 beenVisitedFlag(:) = 1;
                    warning('Clustering aborted since bandwidth far too small! All points in own cluster!');
                    return
                end
            end
            numInitPts_old = numInitPts;
        end
    end
end

