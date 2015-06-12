function close(handles)
    if isfield(handles, 'DH')
        handles.DH.close();
    end