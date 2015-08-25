function [b origInfo] = isZoomed(axesHandle)
    origInfo = getappdata(axesHandle, 'matlab_graphics_resetplotview');
    if isempty(origInfo)
       b = false;
    elseif isequal(get(axesHandle,'XLim'), origInfo.XLim) && ...
           isequal(get(axesHandle,'YLim'), origInfo.YLim) && ...
           isequal(get(axesHandle,'ZLim'), origInfo.ZLim)
       b = false;
    else
       b = true;
    end