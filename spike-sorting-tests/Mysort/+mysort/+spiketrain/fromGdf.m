
function [ST classes ori_classes] = fromGdf(gdf)
    ST = {}; classes = [];
    if isempty(gdf)
        return
    end
    assert(size(gdf,2)==2, 'A gdf has to be a matrix with exactly two colums!');
    ori_classes = unique(gdf(:,1));
    assert(length(ori_classes) < 500, 'This function should not be used for more than 500 neurons!');
    classes = ori_classes;
    zer = find(classes==0,1);
    if ~isempty(zer)
        gdf(:,1) = gdf(:,1) +1;
        classes = classes+1;
    end
    for i=1:length(classes)
        myClass = classes(i);
        ST{i} = gdf(gdf(:,1)==myClass,2); 
    end
end