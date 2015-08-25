function P = spiketrain(gdf, varargin)
P.yoffset = 0;
P.restrictTo = [];
P.colorList = {};
P.setYLims = 0;
P.markersize = 15;
P.plotTicks = 0;
P.linewidth = 1;
P.ah = [];
P.classes = [];
P = mysort.util.parseInputs(P, varargin);

if isempty(gdf)
    return
end

if iscell(gdf) || size(gdf,2) == 1;
    gdf = mysort.spiketrain.toGdf(gdf);
end

if isempty(gdf) % check again
    return
end

if ~isempty(P.restrictTo)
    gdf(~ismember(gdf(:,1), P.restrictTo),:) = [];
end
if isempty(P.classes)
    P.classes = unique(gdf(:,1));
end

if isempty(P.ah)
    P.ah = gca;
end

set(P.ah, 'nextplot', 'add');

for i=1:length(P.classes)
    idx = gdf(:,1) == P.classes(i);
    if sum(idx)==0
        continue
    end
    if isempty(P.colorList)
        myColor = mysort.plot.vectorColor(i);
    else
        myColor = P.colorList{i};
    end
    if P.plotTicks
        st = gdf(idx,2)';
        X = repmat(st,2,1);
        Y = repmat(P.yoffset+i+[0 .8]', 1, length(st));
        line(X,Y, 'parent', P.ah, 'color', myColor, 'linewidth', P.linewidth);
    else
        plot(P.ah, gdf(idx,2), P.yoffset+i, 'o', 'markersize', P.markersize, ...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor', myColor);
    end
end
set(P.ah, 'ytick', (1:length(P.classes))+P.yoffset, 'yticklabel', P.classes)
if P.setYLims
    set(P.ah, 'ylim', [0 length(P.classes)+1]);
end