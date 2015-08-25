function [clustCent,data2cluster,cluster2dataCell] = MeanShiftClusterSpeedBWInc(dataPts,bandWidth,...
            maxBandWidthMultiplicator, nMinPointsInCluster, maxNBandwidthIncreases, bandwidthIncreaseFactor)
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
        if bwIncreaseCounter >= maxNBandwidthIncreases || ~any(removeCluster)
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
            bwIncreaseCounter = bwIncreaseCounter+1;
        end
    end


    %----------------------------------------------------------------------
    function msIteration()
        while numInitPts
            tempInd         = ceil( (numInitPts-1e-6)*rand);        %pick a random seed point
            stInd           = initPtInds(tempInd);                  %use this point as start of mean
            myMean          = dataPts(:,stInd);                           % intilize mean to this points location
            myMembers       = [];                                   % points that will get added to this cluster                          

            sqDistToAll = sum((repmat(myMean,1,numPts) - dataPts).^2);    %dist squared from mean to all points still active

            subidx = find(sqDistToAll < maxBandWidthMultiplicator*bandSq);
            myNumPts = length(subidx);
        %     fprintf('Only %d of %d datapoints in this iteration\n', myNumPts, numPts);

            thisClusterVotes    = zeros(1,myNumPts,'uint16');         %used to resolve conflicts on cluster membership
            firstRound = 1;
            while 1     %loop untill convergence
                if firstRound
                    firstRound = 0;
                    sqDistToAll = sqDistToAll(subidx);
                else
                    sqDistToAll = sum((repmat(myMean,1,myNumPts) - dataPts(:,subidx)).^2);    %dist squared from mean to all points still active
                end
                inInds      = find(sqDistToAll < bandSq);               %points within bandWidth
                thisClusterVotes(inInds) = thisClusterVotes(inInds)+1;  %add a vote for all the in points belonging to this cluster


                myOldMean   = myMean;                                   %save the old mean
                myMean      = mean(dataPts(:,subidx(inInds)),2);                %compute the new mean
                myWeight    = length(inInds);                           % store number of potentially assigned points
                myMembers   = [myMembers inInds];                       %add any point within bandWidth to the cluster
                beenVisitedFlag(subidx(myMembers)) = 1;                         %mark that these points have been visited

                %*** plot stuff ****
        %         if plotFlag
        %             figure(12345),clf,hold on
        %             if numDim == 2
        %                 plot(dataPts(:,subidx)(1,:),dataPts(:,subidx)(2,:),'.')
        %                 plot(dataPts(:,subidx)(1,myMembers),dataPts(:,subidx)(2,myMembers),'ys')
        %                 plot(myMean(1),myMean(2),'go')
        %                 plot(myOldMean(1),myOldMean(2),'rd')
        %                 pause
        %             end
        %         end

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
                        totalWeight = clusterWeights(mergeWith) + myWeight;
                        alpha = clusterWeights(mergeWith)/totalWeight;
                        beta  = myWeight/totalWeight;                        
                        clustCent(:,mergeWith)       = beta*myMean + alpha*clustCent(:,mergeWith);            %record the max as the mean of the two merged (I know biased twoards new ones)
                        clusterVotes(mergeWith,subidx)    = clusterVotes(mergeWith,subidx) + thisClusterVotes;    %add these votes to the merged cluster
                    else    %its a new cluster
                        numClust                    = numClust+1;                   %increment clusters
                        clustCent(:,numClust)       = myMean;                       %record the mean  
                        clusterWeights(1, numClust) = myWeight;
                        clusterVotes(numClust,subidx)    = thisClusterVotes;
                    end

                    break;
                end

            end
            initPtInds      = find(beenVisitedFlag == 0);           %we can initialize with any of the points not yet visited
            numInitPts      = length(initPtInds);                   %number of active points in set
        end
    end
end