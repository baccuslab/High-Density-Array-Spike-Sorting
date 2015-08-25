function handles = SelectUnitIDButton(hObject, eventdata, handles)
    [P PP] = dataviewer.util.getHandleValues(handles);
    id = P.forced_unitID;

    Q = handles.DH.query(['SELECT e.name, b.name, tet.nr, a.id FROM analysis as a '...
                          'JOIN block as b ON (b.id = a.block) '...
                          'JOIN experiment as e ON (e.id = b.expid) '...
                          'JOIN tetrode as tet ON (tet.id = a.tetrode) '...
                          'JOIN unit as u ON (u.analysis = a.id) '...
                          'WHERE u.id = ' num2str(id)]);
            
	if isempty(Q)
        handles.warning_func('UnitID not found!');
        return
    end
    
    figutil.setUserSelection(handles.Experiment, Q{1});
    handles = dataviewer.callback.Experiment(hObject, eventdata, handles);
    
    figutil.setUserSelection(handles.Block, Q{2});
    handles = dataviewer.callback.Block(hObject, eventdata, handles);
    
    figutil.setUserSelection(handles.Tetrode, num2str(Q{3}));
    handles = dataviewer.callback.Tetrode(hObject, eventdata, handles);
    
    set(handles.Analysis, 'Value', find(handles.ids.analysis == Q{4}, 1));
    handles = dataviewer.callback.Analysis(hObject, eventdata, handles);
    
    set(handles.Unit, 'Value', find(handles.ids.units == id, 1));
    handles = dataviewer.callback.Unit(hObject, eventdata, handles);    