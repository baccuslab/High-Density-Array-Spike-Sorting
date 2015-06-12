function TargetSSC = runBOTMOnSpikeSortingContainer(MSDS, sessionIdx, TargetSSC, Tfmax)
    SourceSSC = MSDS.SpikeSortingContainers{MSDS.activeSpikeSortingIdx};
    sI = MSDS.getSessionIterator(sessionIdx);
    cutLeft = SourceSSC.TemplateManager.TemplateList(1).cutLeft;
    TargetSSC.unitNames = SourceSSC.unitNames;
    TargetSSC.TemplateManager = SourceSSC.TemplateManager;    
    
    % There are two options how to prune templates:
    % 1. locally, that means, for every session the templates are pruned
    % for that session only. That means, if we always keep at least 1
    % channel per template, that a session in which a template has only one
    % channel will always keep that channel, even if that channel contains
    % only noise.
    % 2. globally, the templates are pruned only once for the global
    % multi-electrode. That can set a template to zero on complete
    % sessions, even though there might actually be some energy in the
    % template that could be useful for that session
    
    while sI.hasNext()
        DS = sI.next();
        sIdx = sessionIdx(sI.idx);
        ME = DS.MultiElectrode; 
        % SHOULD BE SOLVED already! (Do not try to get sorting from that session. We actually
        % need only the template manager from the multisession sorting
        % object. This will allow to sort also sessions that were not
        % present when the sorting was created.)
        T = SourceSSC.TemplateManager.getWaveforms4MultiElectrode(ME);
        T = mysort.wf.pruneTemplates(T, ...
            'maxChannels', 30,...
            'minChannels', 0, ...
            'absThreshold', 15,...
            'setInvalidChannelsTo', 0);                  
        nC = size(T,2);
        Tf = size(T,1);
        if nargin > 3 && Tf > Tfmax
            Tf = Tfmax;
            T = T(1:Tf, :, :);
        end
            
        noisesd = DS.getNoiseStdFromBufferFile();
        Cest = DS.getCovest();
        botm = mysort.sorters.BOTM(Cest, T,...
            'upsample', 3, 'spikePrior', .01,...
            'max_num_of_channels_per_template', 30);
%         X = DS(1:10000,:);
        gdf = botm.sort(DS);
        gdf(:,2) = gdf(:,2) + cutLeft;
        
        TargetSSC.singleSessionGdfList{sIdx} = gdf;
        TargetSSC.save();
    end  