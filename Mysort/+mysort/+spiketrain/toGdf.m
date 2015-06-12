
function gdf = toGdf(St)
    gdf = [];
    if isempty(St)
        return
    end
%    assert(size(gdf,2)==2, 'A gdf has to be a matrix with exactly two colums!');
    for i=1:length(St)
        gdf = [gdf; [ones(length(St{i}),1)*i mysort.util.toColVec(St{i})]];
    end
    if isempty(gdf)
        return
    end
    gdf = sortrows(gdf,2);
end