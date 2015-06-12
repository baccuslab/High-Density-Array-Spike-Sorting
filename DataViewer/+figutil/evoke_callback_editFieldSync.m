function evoke_callback_editFieldSync(handleName, funh, hObject, eventdata, postCallbackCheck)
    handles = guidata(hObject);
    handles.GUI_CONFIG.(handleName) = funh(get(handles.(handleName), 'string'));
    
    if exist('postCallbackCheck', 'var'),
        handles = postCallbackCheck(handles);
    end
    guidata(hObject, handles);