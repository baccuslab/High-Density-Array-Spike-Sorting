    parfor g=1:nG
       
        if 0
            mysort.plot.waveformsVertical(S)
        end
           
        % Prewhiten
        Spw = mysort.wf.t2v(S)/U;
        % PCA
        F = mysort.util.dimReductionPCA(Spw,...
                P.featureExtraction.nDims, [], 3*1000000);
        % clustering
        maxNBandwidthIncreases = 1;
        bandwidthIncreaseFactor = 1.3;
        maxNClusterSpikes = 200000;
        if size(F,1) <= maxNClusterSpikes
            tic
    %       [clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(X, clustering.bandwidth);
            [clustCent,point2cluster,clustMembsCell] = MeanShiftClusterIncreaseBW(F', P.clustering.meanShiftBandWidth, 0 ,...
                P.clustering.minSpikesPerCluster, maxNBandwidthIncreases, bandwidthIncreaseFactor);
            toc  
            [ids_clu T] = MeanShiftClusterBundleResult(F, clustMembsCell, P.clustering.minSpikesPerCluster);
            cIdx = [];
            T = mysort.util.calculateClassMeans(Spw, ids_clu, 'usemedian');
            if any(unique(ids_clu)==0)
                T = T(2:end,:)';
            end              
        else
            rIdx = randperm(size(F,1));
            cIdx = rIdx(1:min(maxNClusterSpikes, size(F,1)));
            tic
    %       [clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(X, clustering.bandwidth);
            [clustCent,point2cluster,clustMembsCell] = MeanShiftClusterIncreaseBW(F(cIdx,:)', P.clustering.meanShiftBandWidth, 0 ,...
                P.clustering.minSpikesPerCluster, maxNBandwidthIncreases, bandwidthIncreaseFactor);
            toc  
            [ids_clu T] = MeanShiftClusterBundleResult(F(cIdx,:), clustMembsCell, P.clustering.minSpikesPerCluster);
            T = mysort.util.calculateClassMeans(Spw(cIdx,:), ids, 'usemedian');
            if any(unique(ids_clu)==0)
                T = T(2:end,:)';
            end            
        end

        D = bsxfun(@plus, full(dot(Spw',Spw',1))', full(dot(T,T,1))) - full(2*(Spw*T));
        [m ids] = min(D, [], 2);        
        classes = unique(ids);
        if 0 
            [ids(1:10) ids_clu(1:10)]
        end
        CLU{g}.ids = ids;
        CLU{g}.ids_clu = ids_clu;
        CLU{g}.ids_clu_idx = cIdx;
        CLU{g}.spikeSourceFiles = spikeSourceFiles;
        CLU{g}.spikeTimes = spikeTimes;
        CLU{g}.F = F;
        CLU{g}.clusterCenter = clusterCenter;
        CLU{g}.classes = classes;
        CLU{g}.templates = mysort.wf.v2t(mysort.util.calculateClassMeans(mysort.wf.t2v(S), ids, 'usemedian'), nC);
%         clusteredAlignedSpikes = S.spikeAligned.wfs(clustering.clusterIdx,:);
%         clusteredCutSpikes = S.spikeCut.wfs(S.spikeAligned.alignIdx(clustering.clusterIdx),:);
% 
%         clustering.templatesCut = mysort.util.calculateClassMeans(clusteredCutSpikes, clustering.ids, 'usemedian');
%         clustering.templatesAligned = mysort.util.calculateClassMeans(clusteredAlignedSpikes, clustering.ids, 'usemedian');
     
    end  
    disp('Time for big parfor:')
    toc(t1)
disp('Total time:')    
toc(t_glob)

%     for i=1:length(flist)
%         myFile = fullfile(dpath, flist{i});
%     for g=1:nG
%         gdfg = [[CLU(:).ids]
%     
%             
%         gdf = 
%     end


    if 0
        G = ana.mergeLocalSortings(G, meanNoiseStd); 
%         cov(CLU{g}.F(CLU{g}.ids==6,:))        
    end
    if 0

        g = 18;
        figure;
        ah = mysort.plot.subplot2([1 3]);
        mysort.plot.clustering(CLU{g}.F', CLU{g}.ids, [], [], 'axesHandles', ah(2:3));
        mysort.plot.waveformsVertical(CLU{g}.templates, 'IDs', CLU{g}.classes, ...
            'linewidth', 2, 'channelSpaceFactor',1, 'axesHandle', ah(1))
        
        mysort.plot.clustering(CLU{g}.F', CLU{g}.ids_clu, [], []);
    end
    
    if 0 
        %%
        g=1;
        vT = mysort.wf.t2v(CLU{g}.templates);
        nT = size(vT,1);
        CC = zeros(nT,nT);
        for i=1:nT
            for j=i+1:nT
                cc = corrcoef(vT([i j],:)');
                CC(i,j) = cc(2,1);
            end
        end
        figure;
        subplot(1,2,1)
        imagesc(CC)
        colorbar
        subplot(1,2,2)
        imagesc(CC>.9)
        colorbar        
    end