function evoke_callback(fignamespace, callback, hObject, eventdata, postCallbackCheck)
    fname = [fignamespace '.callback.' callback]; 
    callstr = ['handles = ' fname '(hObject, eventdata, guidata(hObject));'];
    eval(callstr);
    if exist('postCallbackCheck', 'var'),
        handles = postCallbackCheck(handles);
    end
    guidata(hObject, handles);