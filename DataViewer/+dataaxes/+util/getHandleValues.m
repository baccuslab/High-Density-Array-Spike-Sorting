function P = getHandleValues(handles)
    P.sliderpos         = get(handles.slider, 'Value');
    P.zoomWindowChecked = get(handles.zoomWindowCheckbox, 'Value');
    P.zoomLevelIdx      = get(handles.zoomPopupmenu, 'Value');
    P.zoomLevelString   = get(handles.zoomPopupmenu, 'String');
    P.zoomLevel         = str2double(P.zoomLevelString(min(P.zoomLevelIdx, length(P.zoomLevelString))));
    P.centerMS          = str2double(get(handles.centralMsEdit, 'String'));
    
    P.dataAxesXLim      = get(handles.dataAxes, 'XLim');
    P.dataAxesYLim      = get(handles.dataAxes, 'YLim');
    