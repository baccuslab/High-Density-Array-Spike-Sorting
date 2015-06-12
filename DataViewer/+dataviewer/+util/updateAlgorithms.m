function handles = updateAlgorithms(handles)
    [P PP] = dataviewer.util.getHandleValues(handles);
    
    % set analysisIDs empty - this is what we want to query
    old_analysisIDs = P.analysisIDs;
    P.analysisIDs = [];
    tidx = []; if ~isempty(PP.trial.val); tidx = str2double(PP.trial.val); end
    analysis = handles.DH.getAnalysis('otherP', P, 'trialIDX', tidx); 
    if isempty(analysis)
        figutil.setListboxString(handles.Unit, '- - -');        
        figutil.setListboxString(handles.Analysis, '- - -');    
        handles.ids.analysis = [];
        return
    end
    names = {}; new_val = [];
    for i=1:size(analysis,1)
        short_descr = analysis{i,2};
        if length(short_descr) > 15
            short_descr = [short_descr(1:12) '...'];
        end
        names{i} = [num2str(analysis{i,1}) ' ' short_descr ' ' analysis{i,4}];
        if old_analysisIDs == analysis{i,1}
            new_val = i;
        end
    end
    figutil.setListboxString(handles.Analysis, names);
    handles.ids.analysis = cell2mat(analysis(:,1));
    if ~isempty(new_val); set(handles.Analysis, 'Value', new_val); end
    [P PP] = dataviewer.util.getHandleValues(handles);
    
    if isempty(old_analysisIDs) || any(P.analysisIDs ~= old_analysisIDs) 
        figutil.setListboxString(handles.Unit, '- - -');
    end
    handles.ids.analysis = cell2mat(analysis(:,1));