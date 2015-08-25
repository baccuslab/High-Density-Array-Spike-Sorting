function handles = PlotXCorrButton(hObject, eventdata, handles)
    [P PP] = dataviewer.util.getHandleValues(handles);
    if isempty(P.trialIDs)
        handles.warning_func('Select Trials first!');
        return
    end
    if isempty(P.unitIDs)
        handles.warning_func('Select Units first!');
        return
    end    
    tGDF = handles.DH.getSpikeTrain('otherP', P);
    spike_trains = mysort.spiketrain.tGdf2cell(tGDF(:,[1 3 5]));
    validTrialIDs = unique(tGDF(:,1));
    trialLength = handles.DH.getTrialLength('trialIDs', validTrialIDs);
    
    srate = handles.DH.getSamplesPerSecond();
    mysort.plot.xcorr(spike_trains, 'srate', srate, 'binSize', .5,...
                      'T', trialLength, 'IDs', PP.unit.names);
    mysort.plot.figureTitle(['Xcorrs Normalized: ' P.figureTitle]);
    handles.new_figure_handles = gcf;  
    mysort.plot.figureName('XCorrs');