function handles = load_3(hObject, eventdata, handles)
    [P PP] = dataviewer.util.getHandleValues(handles);
    
    handles = dataviewer.util.updateTrials(handles);