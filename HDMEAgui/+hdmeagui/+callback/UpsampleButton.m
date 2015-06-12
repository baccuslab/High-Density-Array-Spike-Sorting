function handles = UpsampleButton(hObject, eventdata, handles)
    fprintf('Upsampling Spikes...\n');
    tic
    handles = hdmeagui.data.upsampleSpikes(handles);
    toc
    disp('done');
    