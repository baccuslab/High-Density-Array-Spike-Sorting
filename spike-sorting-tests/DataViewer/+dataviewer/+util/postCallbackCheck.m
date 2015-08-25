function handles = postCallbackCheck(handles)
    P = dataviewer.util.getHandleValues(handles);
    
    % Check the figure handle store for invalid figures
    for i=length(handles.fh_store):-1:1
        if ~ishandle(handles.fh_store(i))
            handles.fh_store(i) = [];
        end
    end
    
    % Treat figure handles, that were created during last callback
    if isfield(handles, 'new_figure_handles') && ...
       ~isempty(handles.new_figure_handles)
        for i=length(handles.new_figure_handles):-1:1
            if ishandle(handles.new_figure_handles(i))
                figutil.attachUserDataToFigure(handles.new_figure_handles(i), P);
            else
                handles.new_figure_handles(i) = [];
            end
        end
        handles.fh_store = [handles.fh_store handles.new_figure_handles];
        handles.new_figure_handles = [];
    end
    
    % Single Trial GROUP
    single_trial_group = {handles.PlotDataButton
                          handles.PlotSpectogram};
    enable = 'on';
    if length(P.trialIDs) ~= 1
        enable = 'off';
    end
    for i=1:size(single_trial_group,1)
        set(single_trial_group{i}, 'Enable', enable);
    end    
    
    % Sorting Button
    if ~isfield(handles, 'dataAxes') || isempty(handles.dataAxes) || ~ishandle(handles.dataAxes)
        set(handles.PlotSortingButton, 'Enable', 'off');
        set(handles.SpikeTrainMode, 'Enable', 'off');      
    else
        set(handles.PlotSortingButton, 'Enable', 'on'); 
        set(handles.SpikeTrainMode, 'Enable', 'on');  
    end
    
    % Trial GROUP
    trial_group = {handles.PlotLoadRewardReaction
                   handles.PlotRates};
    enable = 'on';
    if isempty(P.trialIDs)
        enable = 'off';
    end
    for i=1:size(trial_group,1)
        set(trial_group{i}, 'Enable', enable);
    end

    % Trial - Unit GROUP
    trial_unit_group = {handles.PlotSpikesButton
                        handles.PlotTemplatesButton
                        handles.PlotAmplitudeDriftButton
                        handles.PlotClusteringButton
                        handles.PlotProjectionButton
                        handles.PlotIntraClusterPCAButton
                        handles.PlotISIButton
                        handles.PlotPSTHButton
                        handles.PlotXCorrButton
                        handles.Event
                        handles.Binsize
                        handles.cutleft
                        handles.Tf
                        handles.SpikeMode
                        handles.stacked
                        handles.TemplateModePopup};
    enable = 'on';
    if isempty(P.trialIDs) || isempty(P.unitIDs)
        enable = 'off';
    end
    for i=1:size(trial_unit_group,1)
        set(trial_unit_group{i}, 'Enable', enable);
    end
    
   
    % Covariance GROUP
    covariance_group = {handles.PlotProjectionButton
                        handles.PlotIntraClusterPCAButton
                        handles.PlotClusteringButton};
    if length(P.analysisIDs)~=1
        for i=1:size(covariance_group,1)
            set(covariance_group{i}, 'Enable', 'off');
        end
    end