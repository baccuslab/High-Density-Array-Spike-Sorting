function gdf = gdfMergeClusters(gdf, groups)
    % Merges the clusters in gdf(:,1) accoring to the groups in the cell
    % array groups. The new ID of a group is the index into the groups
    % cell. The old IDs in gdf(:,1) are the values in each group cell
    % element.
    % WARNING! If one or more ID in gdf(:,1) is not in any group it will
    % not be modified. E.g. if gdf(:,1) == 1 is the noise cluster all IDs
    % in groups{1} will be merged into it even if "1" is not in groups{1}
    
    for i=1:length(groups)
        idx = ismember(gdf(:,1), groups{i});
        gdf(idx,1) = -i;
    end
    gdf(:,1) = abs(gdf(:,1));
    