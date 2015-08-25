function handles = trial_rewarded(hObject, eventdata, handles)
    [P PP] = dataviewer.util.getHandleValues(handles);
    
    dataviewer.util.updateTrials(handles);