function handles = PrintSelectedParameters(hObject, eventdata, handles)
    P = dataviewer.util.getHandleValues(handles)
    assignin('base', 'dataViewerP', P);
    
    