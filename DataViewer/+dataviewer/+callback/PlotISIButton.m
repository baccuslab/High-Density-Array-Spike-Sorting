function handles = PlotISIButton(hObject, eventdata, handles)
    P = dataviewer.util.getHandleValues(handles);
    
    GDF = handles.DH.getSpikeTrain('otherP', P);
    gdf = GDF(:,4:5);
    srate = handles.DH.getSamplesPerSecond();
    [spikeTrains IDs] = mysort.spiketrain.fromGdf(gdf);
    mysort.plot.isi(spikeTrains, 'IDs', unique(gdf(:,1)), 'print2ms', 1, 'isHist', 0, 'srate', srate, ...
                    'edges', 0:15:10000, 'figureTitle', P.figureTitle);
	handles.new_figure_handles = gcf;
    mysort.plot.figureName('ISI'); 
    mysort.plot.figureTitle(P.figureTitle);  