function handles = trial_file_error(hObject, eventdata, handles)
    [P PP] = dataviewer.util.getHandleValues(handles);
    
    handles = dataviewer.util.updateTrials(handles);