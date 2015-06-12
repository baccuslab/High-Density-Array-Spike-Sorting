% UPDATE DATA
function DATA = update(handles)
    DATA = handles.DATA;
%     CONFIG = handles.CONFIG;
    GCONFIG = handles.GUI_CONFIG;
    INTER = handles.INTERACTIONS;
    if ~preprocessingOk()
        return
    end

    
    t = INTER.tIdx;
    if INTER.bBrushingChanged
        DATA.bSelectionChanged = 1;
        DATA.templates = hdmeagui.data.templatesBrushingChanged(DATA.templates, t, INTER.bIdx);
    end
    
    if INTER.bBrushingChanged || DATA.bCutSpikesChanged || (isfield(INTER, 'bTemplateChanged') && INTER.bTemplateChanged)
        DATA = hdmeagui.data.cutBrushedSpikes(DATA, INTER, GCONFIG);
        DATA.bCutSpikesChanged = 0;
    end
    DATA = hdmeagui.gui.getSelectedIdx(DATA, INTER, GCONFIG);
    
    if INTER.bBrushingChanged || DATA.bCutSpikesChanged || DATA.bSelectedIdxChanged || DATA.bSpikesAligned
        DATA = hdmeagui.data.calcCurrentTemplate(DATA, INTER, GCONFIG);
    end
    

    %----------------------------------------------------------------------
    function b = preprocessingOk()
        b = 0;
        DATA.bCutSpikesChanged = 0;
        DATA.bSelectionChanged = 0;

        bRawDataPresent   = isfield(DATA, 'X') && isfield(DATA, 'smad') ;
        bCutSpikesPresent = isfield(DATA, 'precutSpikes');
        assert(bRawDataPresent||bCutSpikesPresent, 'Raw data or precut Spikes must be provided!');    

        if bRawDataPresent
            DATA.nC = size(DATA.X,2);
        end
            
        DATA.view_todo.reset = 0;
        DATA.view_todo.updateSelectionPlots = 0;        
        
        bSpikesDetected = isfield(DATA, 'singleChannelDataSorted');
        if ~bSpikesDetected 
            error('not implemented yet');
%             DATA = hdmeagui.data.detectSpikes(DATA, INTER, GCONFIG);
%             DATA = hdmeagui.data.markUselessElectrodes(DATA, INTER, GCONFIG);
            DATA.templates = [];
    %         DATA = hdmeagui.data.cutSpikes(DATA, CONFIG);
            DATA.view_todo.reset = 1;
        end

        bElectrodesMarked = isfield(DATA, 'useElectrodes') && ...
           length(DATA.useElectrodes) == DATA.nC;
        if ~bElectrodesMarked
            DATA = hdmeagui.data.markUselessElectrodes(DATA, INTER, GCONFIG);
    %         DATA = hdmeagui.data.cutSpikes(DATA, GCONFIG);
            DATA.templates = [];
            DATA.view_todo.reset = 1;
        else         
%             % check if spikes are already sufficiently cut
%             bCutSpikesPresent = isfield(DATA, 'cutSpikes') && ...
%                                 isfield(DATA, 'cutSpikesTf') && ...
%                                 isfield(DATA, 'cutSpikesCutleft');
%             bCutSpikesLongEnough = bCutSpikesPresent && ...
%                 DATA.cutSpikesTf >= CONFIG.CutSpikesTf && ...
%                 DATA.cutSpikesCutleft >= CONFIG.CutSpikesCutleft;
%             if ~bCutSpikesLongEnough
%                 DATA = hdmeagui.data.cutBrushedSpikes(DATA, GCONFIG);
%                 DATA.bCutSpikesChanged = 1;
%             end
        end  
        DATA.nUsedElectrodes = sum(DATA.useElectrodes);

        if ~isfield(DATA, 'templates') || isempty(DATA.templates)
            DATA.templates = hdmeagui.data.templatesInit(GCONFIG.CutSpikesTf,...
                 GCONFIG.CutSpikesCutleft, ...            
                 DATA.nUsedElectrodes, size(DATA.singleChannelDataSorted, 1));  
        end
        
        if ~isfield(DATA, 'singleChannelDataModified')
            DATA = hdmeagui.data.calcSpikeAmplitudes(DATA);
        end
         
        b=1;
    end
end

    
