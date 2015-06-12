function handles = Trial(hObject, eventdata, handles)
    % hObject    handle to Trial (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = get(hObject,'String') returns Trial contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from Trial
    [trial, trials, idx, trialID] = figutil.getListBoxValues(handles.Trial, handles.ids.trials);

    if ~isempty(trialID)
        handles = dataviewer.util.updateAlgorithms(handles);
    end
    handles.trialsSelectedBy = 'listbox';