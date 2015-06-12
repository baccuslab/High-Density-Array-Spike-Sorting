function brushedIdx = getBrushedIndex(h)
    if nargin == 0
        h = gca;
    end
    hBrushLine = findall(h, 'tag', 'Brushing');
    brushedData = get(hBrushLine, {'Xdata', 'Ydata'});
    brushedIdx = [];
    if ~isempty(brushedData)
        brushedIdx = ~isnan(brushedData{1})';
    end