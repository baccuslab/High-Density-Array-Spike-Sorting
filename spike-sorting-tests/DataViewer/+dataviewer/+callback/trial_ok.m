function handles = trial_ok(hObject, eventdata, handles)
    [P PP] = dataviewer.util.getHandleValues(handles);
    
    handles = dataviewer.util.updateTrials(handles);