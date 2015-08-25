function handles = Block(hObject, eventdata, handles)
    [block, blocks, idx, blockID] = figutil.getListBoxValues(handles.Block, handles.ids.blocks);
    %dataviewer.util.setListboxString(handles.Tetrode, '- - -');
    %dataviewer.util.setListboxString(handles.Channel, '- - -');
    figutil.setListboxString(handles.Trial, '- - -');
    figutil.setListboxString(handles.Analysis, '- - -');
    figutil.setListboxString(handles.Unit, '- - -');
    handles = dataviewer.util.updateTrials(handles);