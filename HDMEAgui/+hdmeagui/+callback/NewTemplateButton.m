function handles = NewTemplateButton(hObject, eventdata, handles)
    handles.DATA.templates = hdmeagui.data.templatesNew(handles.DATA.templates);
    handles.INTERACTIONS.tIdx = length(handles.DATA.templates.names);
    handles.INTERACTIONS.bTemplateChanged = 1;
    handles.INTERACTIONS.bBrushingChanged = 0;   
    handles.DATA.bCutSpikesChanged = 0;
    handles.DATA.bSelectionChanged = 0;
    handles.DATA.bNeedCutSpikes = 0;

    handles.VIEW = hdmeagui.view.update(handles);   
