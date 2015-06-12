function [gdf_merged, T_merged, localSorting, localSortingID] = computeFinalGDFFromMergedLocalSortings(G)
    gdf_merged = [];
    T_merged = [];
    nextt = 1;
    localSorting = [];
    localSortingID = [];
    for g=1:length(G);
        localValidNeuronsIdx = find(G(g).templates.maxInThisGroup & cellfun(@(x) isempty(x), G(g).templates.masterTemplate));
        if ~isempty(localValidNeuronsIdx)
            localValidNeuronsId = localValidNeuronsIdx-1;
            validGdf = G(g).gdf(ismember(G(g).gdf(:,1), localValidNeuronsId),:); % the -1 is necessary since gdf ids start at zero!
            % set templates without any spikes in gdf as invalid
            validUnits = unique(validGdf(:,1));
            if length(validUnits) < length(localValidNeuronsId)
                warning('There were templates without spikes in the corresponding GDF!! They are removed. This should not happen');
                g
                G(g)
                validUnits
                localValidNeurons
            end
            localValidNeuronsId = intersect(localValidNeuronsId, validUnits);
            localValidNeuronsIdx = localValidNeuronsId+1;
            validGdf(:,1) = validGdf(:,1) + 1000*g;
            gdf_merged = [gdf_merged; validGdf];
            nT = length(localValidNeuronsId);
            T_merged(:,:, nextt:nextt+nT-1) = G(g).templates.wfs(:,:,localValidNeuronsIdx);
            nextt = nextt+nT;
            localSorting = [localSorting; g*ones(nT,1)];
            localSortingID = [localSortingID; localValidNeuronsId(:)];
        end
    end
    units = unique(gdf_merged(:,1));
    nU = length(units);
    assert(length(localSorting) == nU, 'must be identical');
    assert(length(localSortingID) == nU, 'must be identical');
    assert(size(T_merged,3) == nU, 'must be identical');
    assert(nextt-1 == nU, 'must be identical');