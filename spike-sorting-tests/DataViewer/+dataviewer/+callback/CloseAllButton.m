function handles = CloseAllButton(hObject, eventdata, handles)
    % hObject    handle to PlotAmplitudeDriftButton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
   
    for i=length(handles.fh_store):-1:1
        if ishandle(handles.fh_store(i))
            close(handles.fh_store(i));
        end
    end
    handles.fh_store = [];