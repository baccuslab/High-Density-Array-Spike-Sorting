function gdf = mergeUnits(gdf, units)
units = units(:); % make column vector
assert( length(units) > 1, 'More than one units for merging necessary!')

idx = find(ismember( gdf(:,1), units(2:end) )); %find(gdf(:,1) == units(2:end));
gdf(idx ,1) = units(1);
end