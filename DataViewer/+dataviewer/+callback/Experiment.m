function handles = Experiment(hObject, eventdata, handles)
    [experiment, list, idx, experimentID] = figutil.getListBoxValues(handles.Experiment, handles.ids.experiments);

    figutil.setListboxString(handles.Tetrode, '- - -');
    figutil.setListboxString(handles.Block, '- - -');
    figutil.setListboxString(handles.Channel, '- - -');
    figutil.setListboxString(handles.Trial, '- - -');
    figutil.setListboxString(handles.Analysis, '- - -');
    figutil.setListboxString(handles.Unit, '- - -');
    
    if ~isempty(experimentID)
        blocks = handles.DH.getBlocks('experimentIDs', experimentID);
        figutil.setListboxString(handles.Block, blocks(:,2));
        handles.ids.blocks = cell2mat(blocks(:,1));   

        tetrodes = handles.DH.getTetrodes('experimentIDs', experimentID);
        if isempty(tetrodes)
            handles.warning_func(sprintf('This experimentID (%d) has no tetrodes!', experimentID));
            return
        end
        figutil.setListboxString(handles.Tetrode, tetrodes(:,2));
        handles.ids.tetrodes = cell2mat(tetrodes(:,1)); 
    end