function [gdf_merged T_merged localSorting localSortingID G] = processLocalSortings(dpath, runName, elGroupNumbers, elGroupIndices, groupPaths)
    if nargin < 5 groupPaths = []; end

    %% Load GDFs
    
    nGroups = length(elGroupNumbers);
    for i=1:nGroups
        if iscell(groupPaths)
            assert(length(groupPaths) == nGroups, 'Something awefully horrible has happended!');
            subdpath{i} = fullfile(dpath, groupPaths{i});
        else
            subdpath{i} = fullfile(dpath, sprintf('group%04d', i));
            if ~exist(subdpath{i}, 'file')
                subdpath{i} = fullfile(dpath, sprintf('group%03d', i));
            end
        end
    end
    
    G = struct();
    nUnitsPerLocalSorting = [];
    for i=1:nGroups
%         if iscell(groupPaths) & (length(groupPaths) == nGroups)
%             subdpath = fullfile(dpath, groupPaths{i});
%         else
%             subdpath = fullfile(dpath, sprintf('group%03d', i));
%         end
%         S = load(fullfile(subdpath{i}, [runName '.P.mat'])); S=S.S;
        spdet = load(fullfile(subdpath{i}, [runName '.030spikes_det_merged.mat']));
        clust = load(fullfile(subdpath{i}, [runName '.110clusters_meanshift_merged.mat']));
        
        % This line will not work if P.spikeCutting.maxSpikes was set to
        % low in the call for the spike sorting:
        assert(length(clust.clusteringMerged.ids) == length(spdet.spikeDetectionMerged.ts), ...
            ['Not all spikes that were detected were cut (and thus also not matched). This is' ...
             'necessary for the postprocessing. Set P.spikeCutting.maxSpikes higher if you want' ...
             'all spikes to be cut-and-matched.']);
        G(i).gdf = [clust.clusteringMerged.ids spdet.spikeDetectionMerged.ts];
        G(i).units = unique(G(i).gdf(:,1));
        G(i).sortedElectrodes = elGroupIndices{i};
        nUnitsPerLocalSorting(i) = length(G(i).units);
    end
    fprintf('Found Units in local sortings: %d\n', sum(nUnitsPerLocalSorting));

    %% Load originally sorted electrode indices for every group
%     groupfile = fullfile(dpath, [runName '_el_groups.mat']);
%     groups = load(groupfile, 'groups', 'groupsidx', 'nGroupsPerElectrode');
    
    %% Load Templates
    nT = 0;
    for i=1:nGroups
        %subdpath = fullfile(dpath, sprintf('group%03d', i));
        fnames = dir(fullfile(subdpath{i}, [runName '_templates.mat']));
        G(i).templates = load(fullfile(subdpath{i}, fnames.name));
        nT = nT + size(G(i).templates.wfs,3);
    end
    fprintf('Found Templates local sortings: %d\n', nT);

    %% Load Noise
    Ns = cell(nGroups,1);
    for i=1:nGroups
        %subdpath = fullfile(dpath, sprintf('group%03d', i));
        S = load(fullfile(subdpath{i}, [runName '.060cov.mat'])); 
        Ns{i} = S.noise.CestS.CCol;
    end
    meanNoiseStd = mean(sqrt(diag(Ns{1})));

    %% Merge double Templates
    % mergerFile = fullfile(dpath, [runName '_merger.mat']);
    for i=1:length(G)
        assert(length(unique(G(i).gdf(:,1))) == size(G(i).templates.wfs,3), 'must be identical')
    end
    fprintf('Merging double Templates...');
    G = mysort.HDSorting.mergeLocalSortings(G, meanNoiseStd); 
    for i=1:length(G)
        assert(length(unique(G(i).gdf(:,1))) == size(G(i).templates.wfs,3), 'must be identical')
    end    
    g_file = fullfile(dpath, 'G_struct.mat');
    save(g_file, 'G', 'meanNoiseStd');
    
    %% Extract final gdf and templates
    disp('Extracting final gdf and templates...');
    [gdf_merged, T_merged, localSorting, localSortingID] = mysort.HDSorting.computeFinalGDFFromMergedLocalSortings(G);            
    fprintf('Templates after merging: %d\n', size(T_merged,3) ); 
    gdf_merged = sortrows(gdf_merged,2);
    
%     TT = [28 9
%           28 10
%           37 1 
%           11 1];
%     el = G(1).templates.MES.electrodePositions;
%     T = {};
%     for i=1:size(TT,1)
%         T{i} = G(TT(i,1)).templates.wfs(:,:,TT(i,2));
%     end
%     mysort.plot.figure();
%     ah = axes();
%     hold on
%     for i=1:size(TT,1)
%         mysort.plot.templates2D(T{i}, el, 100, 2, 'IDs', i, 'ah', ah)
%     end

%     idx = {};
%     gdfs = {};
%     gdf = [];
%     for i=1:size(TT,1)
%         idx{i} = find(G(TT(i,1)).gdf(:,1) == TT(i,2));
%         gdfs{i} = [i*ones(length(idx{i}),1) G(TT(i,1)).gdf(idx{i},2)];
%         gdf = [gdf; gdfs{i}];
%     end
%     gdf = sortrows(gdf);
%     ovpDist = 1;
%     for i=1:size(TT,1)
%         for j=i+1:size(TT,1);
%             [O nO] = mysort.spiketrain.checkForOverlaps({gdfs{i}(:,2), gdfs{j}(:,2)}, ovpDist);
%             nO
%         end
%     end
%     mysort.plot.xcorr(gdf, 'binSize', .5);
    
