function handles = SelectTrialsetButton(hObject, eventdata, handles)
    % hObject    handle to SelectTrialSetButton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    sel_str = get(handles.Trialset, 'String');
    range = 1:length(handles.ids.trials);
    set(handles.Trial, 'Value', eval(['range(' sel_str ')']));
    handles = dataviewer.util.updateAlgorithms(handles); 
    handles.trialsSelectedBy = 'edit';