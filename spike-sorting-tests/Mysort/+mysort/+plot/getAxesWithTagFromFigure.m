function ah = getAxesWithTagFromFigure(fh, tag)
allAxesInFigure = findall(fh,'type','axes');
ah = allAxesInFigure(arrayfun(@(x) strcmp(get(x, 'Tag'), tag), allAxesInFigure));