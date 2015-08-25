function spiketrain(gdf, varargin)
P.yoffset = 0;
P.restrictTo = [];
P.colorList = {};
P.setYLims = 0;
P.markersize = 15;
P = mysort.util.parseInputs(P, varargin);

if iscell(gdf) || size(gdf,2) == 1;
    gdf = mysort.spiketrain.toGdf(gdf);
end

if ~isempty(P.restrictTo)
    gdf(~ismember(gdf(:,1), P.restrictTo),:) = [];
end
c = unique(gdf(:,1));
hold on
for i=1:length(c)
    idx = gdf(:,1) == c(i);
    if isempty(P.colorList)
        plot(gdf(idx,2), P.yoffset+i, 'o', 'markersize', P.markersize, ...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',mysort.plot.vectorColor(i))
    else
        plot(gdf(idx,2), P.yoffset+i, 'o', 'markersize', P.markersize, ...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',P.colorList{i})
    end
end
set(gca, 'ytick', (1:length(c))+P.yoffset, 'yticklabel', c)
if P.setYLims
    set(gca, 'ylim', [0 length(c)+1]);
end