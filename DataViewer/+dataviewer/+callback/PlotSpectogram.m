function handles = PlotSpectogram(hObject, eventdata, handles)
    P = dataviewer.util.getHandleValues(handles);
    
    [X ] = handles.DH.getData('otherP', P); %, 'downsampleTo', 3000
    srate = handles.DH.getSamplesPerSecond();
    
    % help the fft by truncating x to a length equal to a power of 2
    L = size(X,2);
    L = 2^(nextpow2(L)-1);
    X = X(:, 1:L);
    
    fig = mysort.plot.figure();
    ah = mysort.plot.subplot([size(X,1),1]);
    for i=1:size(X)
        axes(ah(i));
        mysort.util.frequencyAnalysis(X(i,:), srate, 'onesided','.-');
        if i<size(X)
            xlabel('');
        end
        if i>1
            ylabel('');
        end
        title('');
    end
    linkaxes(ah, 'x');
	handles.new_figure_handles = fig;
    mysort.plot.figureName('FrequencyAnalysis'); 
    mysort.plot.figureTitle(P.figureTitle);  