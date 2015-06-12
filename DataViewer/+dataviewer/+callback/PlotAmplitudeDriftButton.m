function handles = PlotAmplitudeDriftButton(hObject, eventdata, handles)
    % hObject    handle to PlotAmplitudeDriftButton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    [P PP] = dataviewer.util.getHandleValues(handles);
    
    [T trialIDs trialidx unitIDs algoIDs channelIDs channelNRs] = ...
        handles.DH.getTemplatesChannel('trialIDs', P.trialIDs, 'unitIDs', P.unitIDs, 'channelIDs', P.channelIDs, 'minmax', true);
    
    mysort.plot.figure()
    chans = unique(channelNRs);
    algos = unique(algoIDs);
    a = mysort.plot.subplot(length(chans));
    for k=1:length(chans)
        idx1 = channelNRs == chans(k);
        axes(a(k));
        hold on
        title(['Channel: ' num2str(chans(k))]);
        for i=1:length(algos)
            idx2 = idx1 & (algoIDs == algos(i));
            plot(trialidx(idx2,1), max(T(idx2,:),[], 2), '-x', 'color', mysort.plot.vectorColor(algos(i)))
            plot(trialidx(idx2,1), min(T(idx2,:),[], 2), '-x', 'color', mysort.plot.vectorColor(algos(i)))
        end
    end
    linkaxes(a, 'xy')
    if min(trialidx)< max(trialidx)
        set(a(1), 'xlim', [min(trialidx) max(trialidx)]);
    else
        set(a(1), 'xlim', [min(trialidx)-1 max(trialidx)+1]);
    end
    set(a(1), 'ylim', [min(T(:)) max(T(:))]);
    xlabel('trial idx')
    ylabel('amplitude')
    mysort.plot.figureTitle(P.figureTitle);
    mysort.plot.figureName('Amplitude Drift'); 
    handles.new_figure_handles = gcf;