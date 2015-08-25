function handles = SaveButton(hObject, eventdata, handles)
    fprintf('Saving...\n');
    
    GUI = handles.GUI;
    CONFIG = handles.CONFIG;
    GUI_CONFIG = handles.GUI_CONFIG;
    version = 1;
    readme = 'This file was created with the meagui for spike sorting. Dont change if you want to reuse it.';
    save('hdmeagui.mat', 'GUI', 'GUI_CONFIG', 'CONFIG', 'version', 'readme', '-v7.3');
    
    disp('Done.');

    