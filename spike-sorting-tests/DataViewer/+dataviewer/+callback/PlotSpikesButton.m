function handles = PlotSpikesButton(hObject, eventdata, handles)
% hObject    handle to PlotSpikesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [P PP] = dataviewer.util.getHandleValues(handles);

    X = handles.DH.getEGDF('otherP', P);
    if isempty(X)
        handles.warning_func('no spikes to plot');
        return
    end
    if P.spikeMode == 1
        mysort.plot.spikes(X(:,3:end), 'classes', X(:,1), 'channelIDs', cell2mat(PP.channel.vals(PP.channel.idx)), 'stacked', P.stacked);
    elseif P.spikeMode == 2
        handles.warning_func('not implemented!');
        return
    end
    handles.new_figure_handles = gcf;
    mysort.plot.figureName('Spikes'); 
    mysort.plot.figureTitle(P.figureTitle);  

