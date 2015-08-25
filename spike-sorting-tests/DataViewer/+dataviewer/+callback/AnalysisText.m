function handles = AnalysisText(hObject, eventdata, handles)
    [ana, anas, idx, anaID] = figutil.getListBoxValues(handles.Analysis, handles.ids.analysis);

    if isempty(anaID)
        return
    end

    [A names]= handles.DH.query(['SELECT * FROM analysis WHERE id IN (' util.dlmstring(anaID) ') ORDER BY date']);
    
%    mysort.plot.printTable(A, 'vlineAfterCols', [1:12], 'colLabel', names);
    mysort.plot.Table('Data', A, 'ColumnName', names);
    handles.new_figure_handles = gcf;