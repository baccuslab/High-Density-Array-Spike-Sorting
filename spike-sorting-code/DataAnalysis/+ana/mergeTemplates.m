function [groups maxT D] = mergeTemplates(T, meanNoiseStd, varargin)
    P.maxRelativeDistance = .35; % in noise std
    P.minCorrelation = .90; % fraction
    P.maxMeanDist = []; %.3; % in noise std
    P = mysort.util.parseInputs(P, varargin, 'error');
    
    [Tf nC nT] = size(T);
    if nT == 1
        groups = {1};
        maxT = T;
        D.MAXDISTS = 0;
        D.DISTS = 0;
        D.Projection = 1;
        D.A = 1;
        D.P = P;
        D.ci = 1;
        return
    end
        
    % only use shift invariant features or use alignment
    MI = squeeze(min(T,[],1));
    MA = squeeze(max(T,[],1));
    M = [MI;MA]/meanNoiseStd;
    
    % remove zero rows
    sM = sum(abs(M), 2);
    M(sM<.1,:) = [];    
    
    % Compute the pairwise similarity measures
    bMerged = true;
    originalIDs = 1:nT;
    newIDs = 1:nT;
    reducedIDs = 1:nT;
    [MAXDISTS DISTS Projection] = computeSimilarities(M);
    MD = MAXDISTS;
    Proj = Projection;
    
    while bMerged
        maxDistCriterion    = MD   < P.maxRelativeDistance;
        projectionCriterion = Proj > P.minCorrelation;          
        A = maxDistCriterion & projectionCriterion & eye(size(MD))==0;
        
        if ~any(A(:))
            bMerged = false;
        else
            MD(~A) = 1;
            [reduced_a reduced_b] = mysort.util.matrixArgMax(-MD);
            a = reducedIDs(reduced_a);
            b = reducedIDs(reduced_b);
            nA = sum(newIDs==a);
            nB = sum(newIDs==b);
            newM = (nA*M(:,reduced_a) + nB*M(:,reduced_b))/(nA+nB);
            assert(nA>0 && nB>0, 'this should not happen!')
            if a<b
                newIDs(newIDs==b) = a;
                M(:,reduced_a) = newM;
                M(:,reduced_b) = [];
                reducedIDs(reduced_b) = [];
            else
                newIDs(newIDs==a) = b;
                M(:,reduced_b) = newM;
                M(:,reduced_a) = [];
                reducedIDs(reduced_a) = [];
            end
            
        end
        [MD DI Proj] = computeSimilarities(M);
    end
    
%     A = sparse(maxDistCriterion & projectionCriterion); 
%     ci = components(A);  
    uGroups = unique(newIDs);
    groups = cell(length(uGroups),1);
    maxT = zeros(1, length(uGroups));
    for i=1:length(uGroups)
        groups{i} = find(newIDs==uGroups(i));
        [m mT] = max(max(abs(mysort.wf.t2v(T(:,:,groups{i}))), [], 2));
        maxT(i) = groups{i}(mT);
    end

    
    %% sort the groups accoring to the id in maxT
    [maxT idx] = sortrows(maxT(:));
    groups_ = groups;
    for i=1:length(idx)
        groups{i} = groups_{idx(i)};
    end
    
    if nargout > 2
        D.MAXDISTS = MAXDISTS;
        D.DISTS = DISTS;
        D.Projection = Projection;
        D.A = A;
        D.P = P;
        D.ci = newIDs;
    end
    
    %----------------------------------------------------------------------
    function [MAXDISTS DISTS Projection] = computeSimilarities(M)
        nT_ = size(M,2);
        MAXDISTS = zeros(nT_);
        DISTS = zeros(nT_);
        Projection = eye(nT_);
        for i=1:nT_
            for j=i+1:nT_
                [MAXDISTS(i,j) idx] = max(abs(M(:,i)-M(:,j)));
                MAXDISTS(i,j) = MAXDISTS(i,j)/max(abs(M(idx,i)),abs(M(idx,j)));
                DISTS(i,j) = sum(abs(M(:,i)-M(:,j)))/size(M,1);
                Projection(i,j) = M(:,i)'*M(:,j)/(norm(M(:,i))*norm(M(:,j)));
                MAXDISTS(j,i) = MAXDISTS(i,j);
                DISTS(j,i) = DISTS(i,j);
                Projection(j,i) = Projection(i,j);
            end
        end
      
    end    
end
    
    
%% Find common groups in adjacency matrix
%     maxDistCriterion = MAXDISTS < .15; % relative to voltage at sample
%     meanDistCriterion = DISTS < .3;   % in noise std
%     projectionCriterion = Projection > 0.92;  % in perc/100
%     
%     A = sparse((maxDistCriterion | meanDistCriterion) & projectionCriterion);% & distCriterion); %distCriterion & corrCriterion);
%     figure;
%     subplot(2,4,1)
%     imagesc(MAXDISTS); colorbar; title('max dist');
%     subplot(2,4,2)
%     imagesc(DISTS); colorbar;  title('mean dist');
%     subplot(2,4,3)
%     imagesc(Projection); colorbar;  title('projection');
%     subplot(2,4,1+4)
%     imagesc(maxDistCriterion); colorbar
%     subplot(2,4,2+4)
%     imagesc(meanDistCriterion); colorbar
%     subplot(2,4,3+4)
%     imagesc(projectionCriterion); colorbar
%     subplot(2,4,4+4)
%     imagesc(A); colorbar    
    
% % %     [xc XD ShiftsXD maxdist maxRelDist] = mysort.wf.tXDiff(resampledTemplates);
% % %     relDistCritetion    = 100*maxRelDist < 30;
% % %     distCriterion = ( maxAbsDistCriterion | relDistCritetion );
% % %     corrCriterion = XC > .9;
% % %     distCriterion = DISTS < 450;
% % %     energyCriterion = E < .1;
    

%     figure;
%     subplot(2,4,1)
%     imagesc(MAXDISTS(sortedgidx,sortedgidx)); colorbar; title('max dist');
%     subplot(2,4,2)
%     imagesc(DISTS(sortedgidx,sortedgidx)); colorbar;  title('mean dist');
%     subplot(2,4,3)
%     imagesc(Projection(sortedgidx,sortedgidx)); colorbar;  title('projection');
%     subplot(2,4,1+4)
%     imagesc(maxDistCriterion(sortedgidx,sortedgidx)); colorbar
%     subplot(2,4,2+4)
%     imagesc(meanDistCriterion(sortedgidx,sortedgidx)); colorbar
%     subplot(2,4,3+4)
%     imagesc(projectionCriterion(sortedgidx,sortedgidx)); colorbar    
%     subplot(2,4,4+4)
%     imagesc(A(sortedgidx,sortedgidx)); colorbar   
    
    
    
%% Or use Clustering    
%         nDim = 14;
%         fetX = mysort.util.dimReductionPCA(M', nDim, [], 3*1000000);
%         bwidth = sqrt(1.5*nDim);
%         [clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(fetX', bwidth);
%         [ids clusterCenter] = MeanShiftClusterBundleResult(fetX, clustMembsCell, 1);
%         uclasses = unique(ids);
%         clustering.templates = mysort.util.calculateClassMeans(VT, ids, 'usemedian');
%         mysort.plot.clustering(fetX, ids);
%         groupIdx = {};
%         groups = {};
%         for i=1:length(uclasses)
%             groupidx = find(ids==uclasses(i));
%             groupIdx{i,1} = groupidx;
%             groups{i,1} = Ugdf(groupidx); 
%             [m maxT] = max(max(abs(mysort.wf.t2v(T(:,:,groupidx))), [], 2));
%             Tmerged(:,:,i) = median(T(:,:,groupidx(maxT)), 3);
%         end
        


%% More fancy features
%     [xc XC ShiftsXC] = mysort.wf.tXCorr(T);