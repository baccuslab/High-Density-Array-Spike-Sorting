function handles = update(handles)
    handles.EVENTS = hdmeagui.gui.update(handles);
    handles.DATA = hdmeagui.data.update(handles);
    handles.VIEW = hdmeagui.view.update(handles);
end