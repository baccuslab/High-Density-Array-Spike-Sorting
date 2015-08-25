function [point relPoint] = getClickPoint(ah)
    point = get(ah, 'CurrentPoint');
    point = point(1,1:2);
    
    if nargout > 1
        xlim = get(ah, 'xlim');
        ylim = get(ah, 'ylim');
        relPoint(1) = (point(1)-xlim(1))/diff(xlim);
        relPoint(2) = (point(2)-ylim(1))/diff(ylim);
    end