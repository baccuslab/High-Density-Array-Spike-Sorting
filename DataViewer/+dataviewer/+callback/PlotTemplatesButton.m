function handles = PlotTemplatesButton(hObject, eventdata, handles)
    P = dataviewer.util.getHandleValues(handles);

    if P.templateMode == 1
        [T, tids, tidx, uids, aids] = handles.DH.getTemplatesConcat('otherP', P);        
    elseif P.templateMode == 2
        [T, aids] = handles.DH.getTemplatesFromCutSpikes('method', 'mean', 'otherP', P);
    elseif P.templateMode == 3
        [T, aids] = handles.DH.getTemplatesFromCutSpikes('method', 'median', 'otherP', P);
    end
    if ~isempty(T) || ~isempty(aids)
        mysort.plot.spikes(T, 'IDs', aids, 'channelIDs', P.channelIDs, 'stacked', P.stacked, 'nC', length(P.channelIDs));
        mysort.plot.figureName('Templates');
        mysort.plot.figureTitle(P.figureTitle);
        handles.new_figure_handles = gcf;        
    else
        warning('no spikes to plot');
    end