function handles = Analysis(hObject, eventdata, handles)
    % hObject    handle to Analysis (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = get(hObject,'String') returns Analysis contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from Analysis
    [ana, anas, idx, anaID] = figutil.getListBoxValues(handles.Analysis, handles.ids.analysis);

    if ~isempty(anaID)
        units = handles.DH.getUnits('analysisIDs', anaID);
        if isempty(units)
            handles.warning_func(sprintf('No units for this Analysis ID (%d)!', anaID));
            return;
        end
        names = {};
        for i=1:size(units,1)
            names{i} = [num2str(units{i,1}) '-' num2str(units{i,2}) '-' units{i,3}];
        end
        figutil.setListboxString(handles.Unit, names);
        handles.ids.units = cell2mat(units(:,1));
        handles.selection_names.unit = cell2mat(units(:,2));
    end