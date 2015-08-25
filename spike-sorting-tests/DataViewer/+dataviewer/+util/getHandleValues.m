function [P PP] = getHandleValues(handles)
    [PP.experiment.val, PP.experiment.vals, PP.experiment.idx, PP.experiment.ids] =  figutil.getListBoxValues(handles.Experiment, handles.ids.experiments);
    [PP.block.val, PP.block.vals, PP.block.idx, PP.block.ids]             = figutil.getListBoxValues(handles.Block, handles.ids.blocks);
    [PP.tetrode.val, PP.tetrode.vals, PP.tetrode.idx, PP.tetrode.ids]     = figutil.getListBoxValues(handles.Tetrode, handles.ids.tetrodes);
    [PP.channel.val, PP.channel.vals, PP.channel.idx, PP.channel.ids]     = figutil.getListBoxValues(handles.Channel, handles.ids.channels);
    [PP.trial.val, PP.trial.vals, PP.trial.idx, PP.trial.ids]             = figutil.getListBoxValues(handles.Trial, handles.ids.trials);
    [PP.unit.val, PP.unit.vals, PP.unit.idx, PP.unit.ids]                 = figutil.getListBoxValues(handles.Unit, handles.ids.units);
    if ~isempty(PP.unit.ids); PP.unit.names = handles.selection_names.unit(PP.unit.idx); else PP.unit.names = []; end
        
    [PP.analysis.val, PP.analysis.vals, PP.analysis.idx, PP.analysis.ids] = figutil.getListBoxValues(handles.Analysis, handles.ids.analysis);
    [PP.event.val, PP.event.vals, PP.event.idx, PP.event.ids] = figutil.getListBoxValues(handles.Event, handles.ids.events);
   
    P.experimentIDs = PP.experiment.ids;
    P.blockIDs      = PP.block.ids;
    P.tetrodeIDs    = PP.tetrode.ids;
    P.channelIDs    = PP.channel.ids;
    P.trialIDs      = PP.trial.ids;
    P.unitIDs       = PP.unit.ids;
    P.analysisIDs   = PP.analysis.ids;
    P.eventIDs      = PP.event.ids;
    
    P.cutleft  = str2double(get(handles.cutleft, 'String'));
    P.Tf       = str2double(get(handles.Tf     , 'String'));
    P.binsize  = str2double(get(handles.Binsize    , 'String'));
    
    P.from     = get(handles.from   , 'String');
    if ~isempty(P.from); P.from = str2double(P.from); else P.from = []; end
    assert(isnumeric(P.from), '"from" must be numeric!');
    if P.from == 1; P.from = []; end
    
    P.to       = get(handles.to     , 'String');
    if isempty(P.to) || strcmp(P.to, 'end'); P.to = []; else P.to = str2double(P.to); end
    assert(isnumeric(P.to), '"to" must be numeric!');
     
    
    P.SpikeTrainMode = figutil.getListBoxValues(handles.SpikeTrainMode);
    P.templateMode = get(handles.TemplateModePopup, 'Value');
    P.spikeMode    = get(handles.SpikeMode, 'Value');
    P.stacked      = get(handles.stacked, 'Value');
    P.trial_ok     = get(handles.trial_ok, 'Value');
    P.trial_file_error = get(handles.trial_file_error, 'Value');
    P.trial_load_error = get(handles.load_error, 'Value');
    P.trial_rewarded = get(handles.trial_rewarded, 'Value');
    P.load = find([get(handles.load_1, 'Value') get(handles.load_2, 'Value') get(handles.load_3, 'Value') get(handles.load_4, 'Value')]);
    P.trialset     = get(handles.Trialset, 'String');
    P.forced_analysisID = str2double(get(handles.ForcedAnalysisID, 'String'));
    P.forced_unitID = str2double(get(handles.ForcedUnitID, 'String'));
    P.trialsSelectedBy = handles.trialsSelectedBy;
    P.figureTitle = dataviewer.util.makeFigureTitle(P, PP);
