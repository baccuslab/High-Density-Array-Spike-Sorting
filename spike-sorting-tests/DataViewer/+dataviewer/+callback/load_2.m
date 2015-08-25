function handles = load_2(hObject, eventdata, handles)
    [P PP] = dataviewer.util.getHandleValues(handles);
    
    handles = dataviewer.util.updateTrials(handles);