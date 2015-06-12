function handles = Tetrode(hObject, eventdata, handles)
    % hObject    handle to Tetrode (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = get(hObject,'String') returns Tetrode contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from Tetrode
    [tetrode, tetrodes, idx, tetrodeIDs] = figutil.getListBoxValues(handles.Tetrode, handles.ids.tetrodes);
    [channel, channels, idx, old_channelID] = figutil.getListBoxValues(handles.Channel, handles.ids.channels);
    if isempty(tetrodeIDs)
        handles.warning_func('No Tetrode selected!');
        return
    end
    channels = handles.DH.getChannels('tetrodeIDs', tetrodeIDs);
    %dataviewer.util.setListboxString(handles.Channel, '- - -');
    figutil.setListboxString(handles.Analysis, '- - -');
    figutil.setListboxString(handles.Unit, '- - -');
    
    if ~isempty(tetrodeIDs)
        handles = dataviewer.util.updateAlgorithms(handles);
        figutil.setListboxString(handles.Channel, channels(:,2));
        if isempty(old_channelID)
            set(handles.Channel, 'Value', [1:size(channels,1)]);
        end
        handles.ids.channels = cell2mat(channels(:,1));
    end