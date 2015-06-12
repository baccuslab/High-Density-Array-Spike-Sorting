function T = templatesBrushingChanged(T, t, bIdx)
    nSp = length(bIdx);
    nSpBrushed = sum(bIdx);
    T.selIdx{t} = true(nSpBrushed,1);
    T.unselIdx{t} = false(nSpBrushed,1);
    if isempty(T.bIdx)
        T.bIdx = false(nSp, t);
    end
    T.bIdx(:, t) = bIdx;
    T.brushed_template(t,:) = -1;
    T.brushed_template_cleaned(t,:) = -1;
    T.selected_template(t,:) = -1;
    T.selected_template_cleaned(t,:) = -1;   