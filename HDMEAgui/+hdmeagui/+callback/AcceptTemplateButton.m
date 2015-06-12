function handles = AcceptTemplateButton(hObject, eventdata, handles)
    fprintf('Accepting Template...\n');
    tic
    handles = hdmeagui.acceptTemplate(handles);
    toc
    disp('done');
    