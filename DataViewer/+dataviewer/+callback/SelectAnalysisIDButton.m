function handles = SelectAnalysisIDButton(hObject, eventdata, handles)
    [P PP] = dataviewer.util.getHandleValues(handles);
    id = P.forced_analysisID;

    Q = handles.DH.query(['SELECT e.name, tet.nr FROM analysis as a '...
                          'JOIN experiment as e ON (a.expid = e.id) '...
                          'JOIN tetrode as tet ON (tet.id = a.tetrode) '...
                          'WHERE a.id = ' num2str(id)]);
            
	if isempty(Q)
        handles.warning_func('AnalysisID not found!');
        return
    end
    
    figutil.setUserSelection(handles.Experiment, Q{1});
    handles = dataviewer.callback.Experiment(hObject, eventdata, handles);
    
    figutil.setUserSelection(handles.Block, '');
    handles = dataviewer.callback.Block(hObject, eventdata, handles);
    
    figutil.setUserSelection(handles.Tetrode, num2str(Q{2}));
    handles = dataviewer.callback.Tetrode(hObject, eventdata, handles);
    
    set(handles.Analysis, 'Value', find(handles.ids.analysis == id, 1));
    handles = dataviewer.callback.Analysis(hObject, eventdata, handles);