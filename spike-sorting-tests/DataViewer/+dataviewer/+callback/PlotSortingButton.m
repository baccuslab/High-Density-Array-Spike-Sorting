function handles = PlotSortingButton(hObject, eventdata, handles)
    P = dataviewer.util.getHandleValues(handles);
    if length(P.trialIDs) ~= 1; warning('Only a single trial can be plotted at the moment!'); return; end
  
    if ~isfield(handles, 'dataAxes') || ~ishandle(handles.dataAxes)
       handles.warning_func('Create a data figure first!')
       return
    end
    
    gdf = handles.DH.getGDF('otherP', P);
    if isempty(gdf)
        handles.warning_func('No Spikes To Plot!')
        return
    end   
    
    if strcmp(P.SpikeTrainMode, 'Templates')
        %P.channelIDs = []; P.channelNRs = [];
        [T, trialidx, trialIDs, unitIDs, algoID, cutleft] = handles.DH.getTemplatesConcat('otherP', P);
        if isempty(T) || size(T,1) < length(P.unitIDs)
            handles.warning_func('At least one selected Units has no Template!!!');
            return
        end
        nC = 4;
  
        gdf(:,2) = gdf(:,2)-(cutleft(1)); 
      
        %mysort.plot.sorting('axesHandle', gca, 'spacer', dataAxesSpacer, 'sampleOffset', P.from, 'srate', handles.DH.getSamplesPerSecond());
        %axes(handles.dataAxes);
        T = mysort.wf.v2t(T, nC);
        T = T(:, size(T,2):-1:1, :);
        mysort.plot.templateSpikeTrain(T, gdf, 'channelSpacer', handles.dataAxesSpacer, ...
            'AxesHandle', handles.dataAxes, 'mode', 'normal', ...
            'timeMultiplicator', handles.DH.getSrate()*1000, ...
            'T_gdf_idx2id', algoID);
    elseif strcmp(P.SpikeTrainMode, 'Lines')
        axes(handles.dataAxes);
        classes = unique(gdf(:,1));
        for i=1:length(classes)
            mysort.plot.verticalLines(gdf(gdf(:,1)==classes(i),2)/(handles.DH.getSamplesPerSecond()/1000), [], 'color', mysort.plot.vectorColor(classes(i)));
        end
    elseif strcmp(P.SpikeTrainMode, 'Dots')
        handles.warning_func('not implemented');
    end