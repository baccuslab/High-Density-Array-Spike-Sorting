function handles = AlignButton(hObject, eventdata, handles)
    fprintf('Aligning Spikes...\n');
    tic
    handles = hdmeagui.GUI.alignSpikes(handles);
    toc
    disp('done');