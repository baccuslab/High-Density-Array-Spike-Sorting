function localSortings = createUMS2000Structure(dpath, runName, hdmea, sessionIdx,...
                                        sessionLengths, el_groups, ME,...
                                        T, localIdx, localID)
    
    samplesPerSecond = hdmea.getSamplesPerSecond();
    nT = size(T,3);
    nLS = length(el_groups);
    localSortings = cell(1,nLS);
    if isa(hdmea, 'mysort.ds.MultiSessionInterface')
        all_session_fileNames = hdmea.getSessionFilenames();
        configME = hdmea.getMergedMultiElectrode4Sessions(sessionIdx);
    else
        if ismember(hdmea, 'fname')
            all_session_fileNames = hdmea.fname;
        else
            all_session_fileNames = 'bla';
        end
        configME = hdmea.MultiElectrode;
    end
    
    usesSessionFileNames = all_session_fileNames(sessionIdx);
    
    light_ts = getLightTS();
    
    xy = configME.getElectrodePositions();
    configMEIdxForGlobalME = find(ismember(ME.electrodeNumbers, configME.electrodeNumbers));
    cutStart = 4;
    TfMichele = 30;
    cutStart = 1;
    TfMichele = size(T,1);
    T = T(cutStart:cutStart+TfMichele-1,configMEIdxForGlobalME,:);
        
    for myLocalSortingIdx=1:nLS
        n = struct();
        myLocalNeuronsIdx = find(localIdx==myLocalSortingIdx);
        elNumbers = el_groups{myLocalSortingIdx};
        [b elIndices] = ismember(elNumbers, configME.electrodeNumbers);
        nC = length(elNumbers);
        
        myLocalAcceptedClusterIDs = localID(myLocalNeuronsIdx)+1; % shift by 1 for the 0 noise cluster
        subdpath = fullfile(dpath, sprintf('group%03d', myLocalSortingIdx));
        P = load(fullfile(subdpath, [runName '.P.mat']));
%         DET_raw = load(fullfile(subdpath, [runName '.020spikes_det.mat']));
        DET = load(fullfile(subdpath, [runName '.030spikes_det_merged.mat']));
        SP = load(fullfile(subdpath, [runName '.040spikes_cut.mat']));
        N  = load(fullfile(subdpath, [runName '.060cov.mat']));
%             AL = load(fullfile(subdpath, [runName '.050spikes_aligned.mat']));
%             CL = load(fullfile(subdpath, [runName '.090clusters_meanshift.mat']));
        CLTM = load(fullfile(subdpath, [runName '.100botm_matching.mat']));
        CLM = load(fullfile(subdpath, [runName '.110clusters_meanshift_merged.mat']));        
%             spikesUsedForClustering = SP.spikeCut.cutIdx(AL.spikeAligned.alignIdx(CL.clustering.clusterIdx));
        myTs_unwrapped = DET.spikeDetectionMerged.ts(SP.spikeCut.cutIdx)/samplesPerSecond;
        myTs = mysort.spiketrain.splitGdf(myTs_unwrapped, [0 cumsum(sessionLengths/samplesPerSecond)]);
        myTrials = cell2mat(arrayfun(@(x,y) ones(x,1)*y,cellfun(@length, myTs), 1:length(myTs),'UniformOutput', false)');
        myTs = cell2mat(myTs');
        localSortingAssignments = CLM.clusteringMerged.ids+1; % shift because of the zero cluster
        myClusterIDs = unique(localSortingAssignments);
        
        n.x = xy(:,1)';
        n.y = xy(:,2)';
        % Store the electrode numbers of the whole electrode configuration
        n.el_idx = configME.electrodeNumbers;
        % Now store the index into the electrode number list for the electrodes 
        % that were used in the local sorting of that neuron were used 
        n.clustered_source = elIndices;
        
        n.trace_name = usesSessionFileNames;
        n.clus_of_interest = myLocalAcceptedClusterIDs;
        
        n.light_ts = light_ts;

        n.spikes.params = struct();
        n.spikes.clus_of_interest = n.clus_of_interest;
        for i=1:length(myLocalNeuronsIdx)
            n.spikes.template_of_interest{n.clus_of_interest(i)} = T(:,:,myLocalNeuronsIdx(i));
            n.pktk(:,i) = max(T(:,:,myLocalNeuronsIdx(i)),[],1) - min(T(:,:,myLocalNeuronsIdx(i)),[],1);
        end
        n.spikes.x = n.x;
        n.spikes.y = n.y;
        if isempty(CLTM.clusteringMatched.spikeCutAligned)
            twfs = mysort.wf.v2t(SP.spikeCut.wfs, nC);
        else
            twfs = mysort.wf.v2t(CLTM.clusteringMatched.spikeCutAligned, nC);
        end
        twfs = twfs(cutStart:min(end, cutStart+TfMichele-1),:,:);
        n.spikes.waveforms = permute(twfs,[3 1 2]);
        n.spikes.spiketimes = myTs';
        n.spikes.unwrapped_times = myTs_unwrapped';
        n.spikes.assigns = localSortingAssignments'; 
        n.spikes.trials = myTrials';
        
        % Set the spike labels for the cluster of interest,i.e., the
        % accepted neurons to accepted.
        n.spikes.labels = [(1:max(myClusterIDs))' ones(max(myClusterIDs),1)];
        n.spikes.labels(ismember(myClusterIDs, myLocalAcceptedClusterIDs),2) = 2; % good units
        n.spikes.labels(1,2) = 4; % garbage unit that was "0" before
        
        % Info Struct
        n.spikes.info = struct();
        n.spikes.info.detect.stds = sqrt(diag(N.noise.CestS.CCol(1:nC,1:nC)));
        n.spikes.info.detect.dur = sessionLengths/samplesPerSecond;
        n.spikes.info.detect.thresh = (-3.5*sqrt(diag(N.noise.CestS.CCol(1:nC,1:nC))))';% P.S.P.spikeDetection.thr;
        n.spikes.info.pca = [];
        n.spikes.info.align = [];
        n.spikes.info.kmeans.colors = getCMap();
        n.spikes.info.kmeans.assigns = localSortingAssignments';
        n.spikes.info.kmeans.num_clusters = max(myClusterIDs);
        n.spikes.info.interface_energy = [];
        n.spikes.info.tree = [];
        % Params Struct
        n.spikes.params.Fs = samplesPerSecond;
        % the following values come from one of Micheles sortings
        n.spikes.params.detect_method = 'auto';
        n.spikes.params.thresh = 4;
        n.spikes.params.window_size = 1.5;
        n.spikes.params.shadow = .75;
        n.spikes.params.cross_time = .6;
        n.spikes.params.refractory_period = 1.5;
        n.spikes.params.max_jitter = .6;
        n.spikes.params.agg_cutoff = .05;
        n.spikes.params.kmeans_clustersize = 500;
        % Display struct
        n.spikes.params.display.aspect_ratio = .667; 
        n.spikes.params.display.default_waveformmode = 2; 
        n.spikes.params.display.time_scalebar = 1; 
        n.spikes.params.display.cmap = getCMap(); 
        n.spikes.params.display.xchoice = 'PC'; 
        n.spikes.params.display.xparam = 1; 
        n.spikes.params.display.ychoice = 'PC'; 
        n.spikes.params.display.yparam = 2; 
        n.spikes.params.display.show_outliers = 1; 
        n.spikes.params.display.show_isi = 1; 
        n.spikes.params.display.max_autocorr_to_display = .1; 
        n.spikes.params.display.max_isi_to_display = .025; 
        n.spikes.params.display.correlations_bin_size = 2; 
        n.spikes.params.display.isi_bin_size = .1; 
        n.spikes.params.display.default_xcorr_mode = 1; 
        n.spikes.params.display.trial_spacing = .5; 
        n.spikes.params.display.stability_bin_size = .5; 
        n.spikes.params.display.max_scatter = 1000; 
        n.spikes.params.display.default_outlier_method = 1; 
        n.spikes.params.display.label_categories = {'in process','good unit','multi-unit','garbage','needs outlier removal';}; 
        n.spikes.params.display.label_colors = [0.700000000000000,0.700000000000000,0.700000000000000;0.300000000000000,0.800000000000000,0.300000000000000;0.300000000000000,0.300000000000000,0.800000000000000;0.800000000000000,0.300000000000000,0.300000000000000;0.700000000000000,0.700000000000000,0.300000000000000;];
        n.spikes.params.display.default_figure_size = [0.050000000000000,0.100000000000000,0.900000000000000,0.800000000000000]; 
        n.spikes.params.display.figure_font_size = 8; 
        n.spikes.params.display.initial_split_figure_panels = 4; 
        n.spikes.params.display.merge_fig_color = [0.700000000000000,0.800000000000000,0.700000000000000]; 
        n.spikes.params.display.split_fig_color = [0.800000000000000,0.700000000000000,0.700000000000000]; 
        n.spikes.params.display.outlier_fig_color = [0.700000000000000,0.700000000000000,0.800000000000000]; 
        n.spikes.params.display.margin = 140; 
        n.spikes.params.display.outer_margin = 80; 
        n.spikes.params.display.width =140;
        
        localSortings{myLocalSortingIdx} = n;
        fprintf('%d%% ', round(100*myLocalSortingIdx/nLS));
    end
    fprintf('\n');
    
    %----------------------------------------------------------------------
    function lts = getLightTS()
        lts = {};        
    end
    %----------------------------------------------------------------------
    function cm = getCMap()
        cm = [0.0416666666666667,0,0;0.0833333333333333,0,0;0.125000000000000,0,0;0.166666666666667,0,0;0.208333333333333,0,0;0.250000000000000,0,0;0.291666666666667,0,0;0.333333333333333,0,0;0.375000000000000,0,0;0.416666666666667,0,0;0.458333333333333,0,0;0.500000000000000,0,0;0.541666666666667,0,0;0.583333333333333,0,0;0.625000000000000,0,0;0.666666666666667,0,0;0.708333333333333,0,0;0.750000000000000,0,0;0.791666666666667,0,0;0.833333333333333,0,0;0.875000000000000,0,0;0.916666666666667,0,0;0.958333333333333,0,0;1,0,0;1,0.0416666666666667,0;1,0.0833333333333333,0;1,0.125000000000000,0;1,0.166666666666667,0;1,0.208333333333333,0;1,0.250000000000000,0;1,0.291666666666667,0;1,0.333333333333333,0;1,0.375000000000000,0;1,0.416666666666667,0;1,0.458333333333333,0;1,0.500000000000000,0;1,0.541666666666667,0;1,0.583333333333333,0;1,0.625000000000000,0;1,0.666666666666667,0;1,0.708333333333333,0;1,0.750000000000000,0;1,0.791666666666667,0;1,0.833333333333333,0;1,0.875000000000000,0;1,0.916666666666667,0;1,0.958333333333333,0;1,1,0;1,1,0.0625000000000000;1,1,0.125000000000000;1,1,0.187500000000000;1,1,0.250000000000000;1,1,0.312500000000000;1,1,0.375000000000000;1,1,0.437500000000000;1,1,0.500000000000000;1,1,0.562500000000000;1,1,0.625000000000000;1,1,0.687500000000000;1,1,0.750000000000000;1,1,0.812500000000000;1,1,0.875000000000000;1,1,0.937500000000000;1,1,1;];
    end
end