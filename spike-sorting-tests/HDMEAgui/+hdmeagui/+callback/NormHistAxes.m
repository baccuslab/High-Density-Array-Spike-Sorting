function handles = NormHistAxes(hObject, eventdata, handles)
    point = get(gca, 'CurrentPoint');
    handles.GUI.normHistThresh = point(1);
    handles.GUI.normHistThreshWasChanged = 1;
    if ishandle(handles.GUI.normHistPlotThresholdHandle)
      delete(handles.GUI.normHistPlotThresholdHandle);
    end
    handles.GUI.normHistPlotThresholdHandle = ...
        mysort.plot.verticalLines(point(1), [], 'color', handles.CONFIG.color_normhist_new_thresh);
    handles = handles.GUI.update(handles);