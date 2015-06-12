function handles = updateTrials(handles)
    [P PP] = dataviewer.util.getHandleValues(handles);
    
    if ~isempty(P.blockIDs)
        if ~P.trial_ok
            P.trial_ok = [];
        end
        if P.trial_rewarded == 1
            rew = [];
        elseif P.trial_rewarded == 2
            rew = true;
        else
            rew = false;
        end
        
        trials = handles.DH.getTrials('blockIDs', P.blockIDs, 'ok', P.trial_ok, ...
            'rewarded', rew, 'ignoreFileErrorTrials', P.trial_file_error, ...
            'ignoreLoadErrorTrials', P.trial_load_error, 'load', P.load);
        if ~isempty(trials)
            figutil.setListboxString(handles.Trial, trials(:,2));
            handles.ids.trials = cell2mat(trials(:,1));
        else
            figutil.setListboxString(handles.Trial, '- - -');
            handles.ids.trials = [];            
        end
        
        handles = dataviewer.util.updateAlgorithms(handles);
    end