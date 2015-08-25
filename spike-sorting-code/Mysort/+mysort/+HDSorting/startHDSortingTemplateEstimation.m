function startHDSortingTemplateEstimation(outPath, dpath, runName, Tf, cutleft, DSFull, elGroupIndices, MES)  %elGroupIndices is necessary for saving!
    % This function is only necessary since matlab cannot save in a par for
    % loop if the save is called in the first layer of the body of the loop
    matchingFile  = fullfile(outPath, dpath, [runName '.100botm_matching.mat']);
    sortingFile   = fullfile(outPath, dpath, [runName '.110clusters_meanshift_merged.mat']);
    templateFile  = fullfile(outPath, dpath, [runName '_templates.mat']);        

    if exist(templateFile, 'file')
        disp('Templates already computed, skipping');
        return
    end
    
    S = load(sortingFile);
    M = load(matchingFile);
    units = unique(S.clusteringMerged.ids);
    wfs = zeros(Tf, size(DSFull,2), length(units));
    nSourceSpikesPerTemplateAndChannel = [];
    
    % Select from all neurons maximally 300 spikes and cut them all at once
    % %%
    % WARNING This gdf has its two columns switched !!!
    % %% 
    igdf = zeros(length(units)*300, 2);
    lastIdx = 0;
    for uidx = 1:length(units)
        tsIdx = find(S.clusteringMerged.ids == units(uidx));
        tsIdx = tsIdx(1:min(300, end));
        nS = length(tsIdx);
        igdf(lastIdx+1:lastIdx+nS,:) = [round(M.clusteringMatched.ts(tsIdx)) ones(nS,1)*uidx];
        lastIdx = lastIdx + nS;
    end
    igdf(lastIdx+1:end,:) = [];
    igdf = sortrows(igdf);
    
    if ~isempty(igdf)
        fprintf('Cutting Spikes on all channels to estimate full templates...\n');
        t = tic;
        wfs_ = DSFull.getWaveform(igdf(:,1), cutleft, Tf);
        t = toc(t);
        fprintf('Cutting Spikes done. (%.2fs)\n', t);
    
        for uidx = 1:length(units)
            myidx = igdf(:,2) == uidx;
            wfs(:,:,uidx) = mysort.wf.v2t(median(wfs_(myidx,:),1), size(DSFull,2));
            nSourceSpikesPerTemplateAndChannel(uidx, 1:size(DSFull,2)) = size(wfs_,1); % necessary for saving!
        end
    else
        wfs = [];
        nSourceSpikesPerTemplateAndChannel = [];
    end
    
    save(templateFile, 'wfs', 'cutleft', 'MES', 'elGroupIndices', 'nSourceSpikesPerTemplateAndChannel'); 
    