function title = makeFigureTitle(P, PP)
    expstr = ['E:' PP.experiment.val{1}];
    
    blockstr = [' B:' PP.block.val{1}];
    for i=2:length(PP.block.val)
        blockstr = [blockstr ' ' PP.block.val{i}];
    end
    
    tstr = '  Trial:';
    trialsetstr = '';
    if strcmp(P.trialsSelectedBy, 'listbox')
        sel_str = dataviewer.util.getStringFromSelection(PP.trial.val);
        trialsetstr = [tstr sel_str];
    elseif strcmp(P.trialsSelectedBy, 'edit')
        trialsetstr = [tstr P.trialset];
    end
    
    tetstr = ['  Tet:' dataviewer.util.getStringFromSelection(PP.tetrode.val)];
    
    anastr = [];
    if ~isempty(P.analysisIDs)
        anastr = ['  A:' dataviewer.util.getStringFromSelection(P.analysisIDs)];
    end
    
    unitstr = [];
    if ~isempty(P.unitIDs)
        unitstr = ['  U:' dataviewer.util.getStringFromSelection(P.unitIDs)]; %PP.unit.names
    end
    title = [expstr blockstr trialsetstr tetstr anastr unitstr];
    