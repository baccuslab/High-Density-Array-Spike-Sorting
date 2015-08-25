function openSplitMergeTool(ename, ntkfilelist, spikeSortingRunName, elGroups)
    spsort_path = '/links/groups/hima/recordings/HiDens/SpikeSorting/';
    egroup = 'Roska';

    if ~iscell(ntkfilelist)
        ntkfilelist = {ntkfilelist};
    end
    
    % set some paths
    spsort_out_path = fullfile(spsort_path, egroup, ename);

    for k=1:length(ntkfilelist)
        [PATH,NAME,EXT] = fileparts(ntkfilelist{k});
        spath = fullfile(spsort_out_path, NAME, 'sortings', spikeSortingRunName); 
        export_file = fullfile(spath, [spikeSortingRunName 'Export4UMS2000']);
        R = load(export_file);    
        for g=1:length(elGroups)
            splitmerge_tool(R.rgcs{elGroups(g)}.spikes);
        end
    end





%% This is old code to check if the merging worked properly
% for k=2:length(files_to_open_spit_merge_tool)
%     [PATH,NAME,EXT] = fileparts(ntkfilelist{files_to_open_spit_merge_tool(k)});
%     ntkfile = [NAME EXT];
%     hd5file = fullfile(spsort_out_path, [NAME '.h5']);
%     gdffile = fullfile(spsort_out_path, [NAME '.gdf']);
%     spath = fullfile(spsort_out_path, NAME, 'sortings');    
%     grouppath = fullfile(spath, sprintf('group%03d', selected_electrode_group));
%     export_file = fullfile(spath, [spikeSortingRunName 'Export4UMS2000']);
%     R{k} = load(export_file);
%     splitmerge_tool(R{k}.rgcs{selected_electrode_group}.spikes);
%     
%     S = struct();
%     load(fullfile(grouppath, [spikeSortingRunName '.P.mat']));
%     P = S.P;    
%     a = load(fullfile(grouppath, [spikeSortingRunName '.040spikes_cut.mat']));
%     S.spikeCut = a.spikeCut;
%     a = load(fullfile(grouppath, [spikeSortingRunName '.050spikes_aligned.mat']));
%     S.spikeAligned = a.spikeAligned;
%     a = load(fullfile(grouppath, [spikeSortingRunName '.060cov.mat']));
%     S.noise = a.noise;
%     a = load(fullfile(grouppath, [spikeSortingRunName '.090clusters_meanshift.mat']));
%     S.clustering = a.clustering;
%     a = load(fullfile(grouppath, [spikeSortingRunName '.100botm_matching.mat']));
%     S.clusteringMatched = a.clusteringMatched;
%     a = load(fullfile(grouppath, [spikeSortingRunName '.110clusters_meanshift_merged.mat']));
%     S.clusteringMerged = a.clusteringMerged;
%     
%     clear a;
%     
%     mysort.plot.waveforms(S.spikeAligned.wfs, 'nC', nC, 'IDs', S.clustering.ids, 'stacked', 0)
%     mysort.plot.waveforms(S.spikeAligned.wfs, 'nC', nC, 'IDs', S.clusteringMatched.ids, 'stacked', 0)
%     mysort.plot.waveforms(S.spikeAligned.wfs, 'nC', nC, 'IDs', S.clusteringMerged.ids, 'stacked', 0)
% 
%     S.clusteringMatched.template_ids = unique(S.clusteringMatched.ids);
%     % in S.clusteringMatched.ids the ids are going from 0 to nT with
%     % 0 being not matched and noise, thus there is one "real" template less
%     % than units (if any spike has an ID of 0 !!)
%     realTemplateIdx = find(S.clusteringMatched.template_ids>0);
%     nT = length(realTemplateIdx);
%     if nT == 0
%         % we have only the noise cluster
%         D = [];
%         maxT = [];
%         groups = {};
%     else
%         % resample only templates that do not have id==0
%         resampledTemplates = mysort.util.resampleTensor(mysort.wf.v2t(...
%             S.clusteringMatched.templates(realTemplateIdx,:), nC), 3, 1);
%         mysort.plot.waveforms(resampledTemplates, 'IDs', S.clusteringMatched.template_ids(realTemplateIdx), 'stacked', 0);
%         
%         % the indices in the groups are indices into
%         % realTemplateIdx which is index into S.clusteringMatched.template_ids
%         [groups maxT D] = ana.mergeTemplates(resampledTemplates, S.noise.meanNoiseStd,...
%             'maxRelativeDistance', P.mergeTemplates.ifMaxRelDistSmallerPercent/100,...
%             'minCorrelation', P.mergeTemplates.atCorrelation);
%     end
%     clusteringMerged.D = D;
%     clusteringMerged.maxT = maxT;
%     clusteringMerged.groups = groups;            
%     clusteringMerged.templates = zeros(length(groups), P.spikeCutting.Tf*nC);
%     clusteringMerged.ids = zeros(length(S.clusteringMatched.ids),1);
%     
%     % we need to merge now the indices in ids accoring to the groups in
%     % groups. But groups points into realTemplateIdx !!
%     for i=1:length(groups)
%         group_ids = S.clusteringMatched.template_ids(realTemplateIdx(clusteringMerged.groups{i}));
%         idx = ismember(S.clusteringMatched.ids, group_ids);
%         clusteringMerged.ids(idx) = i;
%         clusteringMerged.templates(i,:) = median(S.spikeCut.wfs(idx,:),1);
%     end 
%     mysort.plot.waveforms(S.spikeAligned.wfs, 'nC', nC, 'IDs', clusteringMerged.ids, 'stacked', 0)
% end
% 

