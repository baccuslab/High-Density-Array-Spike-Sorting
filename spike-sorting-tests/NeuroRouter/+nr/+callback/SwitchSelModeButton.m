function gui = SwitchSelModeButton(hObject, eventdata, gui)
    if strcmp(gui.selectionMode, 'add')
        gui.selectionMode = 'remove';
        set(gui.SelModeButton, 'String', 'Mode: Remove');
    else
        gui.selectionMode = 'add';
        set(gui.SelModeButton, 'String', 'Mode: Add');
    end