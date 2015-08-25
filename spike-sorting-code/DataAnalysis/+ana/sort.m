function [S P] = sort(DS, dpath, name, varargin)
    warning('This function is depricated, use mysort.sorters.sort() instead');
    error('Depricated!')
    % Artefact Detection
    P.artefactDetection.use = 0;
    P.artefactDetection.width = 350;
    P.artefactDetection.threshold = size(DS,2)*41; % should be sampling frequency dependent!
    
    % Spike Detection
    P.spikeDetection.method = '-';
    P.spikeDetection.thr = 3.5;
    P.spikeDetection.minDist = 25;
    P.spikeDetection.maxDataLength = [];
    P.spikeDetection.mergeSpikesMaxDist = 18; % do not go below 13 for HDMEA!
    P.spikeDetection.removeEventsWithAbsAmplitudeLargerThan = 15000;
    P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude = [0 0 0]; % should remove atrifacts that are oscillations (successive detections of peaks with alternating signs). First parameter is N, second the minHeight, third the maxDist within a group
    
    % Spike cutting
    P.spikeCutting.maxSpikes = 100000;
    P.spikeCutting.Tf = 55;
    P.spikeCutting.cutLeft = 10;

    % Spike Alignment
    P.spikeAlignment.method = 'onUpsampledMean';
    P.spikeAlignment.maxSpikes = 60000;
    P.spikeAlignment.Tf = 25;
    P.spikeAlignment.cutLeft = 6;
    P.spikeAlignment.initAlignment = '-';
    P.spikeAlignment.maxIdx = P.spikeAlignment.cutLeft + 1;
    
    P.spikeAlignment.maxIterations = 30;
    
    % Noise estimation    
    P.noiseEstimation.minLength = 2000000;
    P.noiseEstimation.minDistFromSpikes = 60;
    
    % Feature extraction
    P.featureExtraction.Tf = 15;
    P.featureExtraction.cutLeft = 3;
    P.featureExtraction.nDims = 6;
    
    % Clustering
    P.clustering.maxSpikes = 30000;
    P.clustering.meanShiftBandWidth = []; % default is sqrt(1.1*nDims)
%     P.clustering.bandwidthMultMaxDistSpeedMS = 5;
    P.clustering.minSpikesPerCluster = 10;
    
    % Template Matching on Cut Spikes
    P.templateMatchingCut.prior = .0001;
    P.templateMatchingCut.residualPeakSmallerThanStdNoise = [];
    
    % Merge Templates after template matching
    P.mergeTemplates.merge = 1;
	P.mergeTemplates.upsampleFactor = 3;
    P.mergeTemplates.atCorrelation = .96;
%     P.mergeTemplates.ifMeanDistSmaller = .1;    % in units std of noise
    P.mergeTemplates.ifMaxDistSmaller = 2.5;      % in units std of noise
    P.mergeTemplates.ifMaxRelDistSmallerPercent = 25;

    % Run BOTM Template Matching on whole data
    P.botm.run = 0;
    P.botm.Tf = 55;
    P.botm.cutLeft = 10;
    P.botm.prior = .0001;

    P = mysort.util.parseInputs(P, varargin, 'error');

    % Checks
    assert(P.botm.run==0 ||P.botm.cutLeft <= P.spikeCutting.cutLeft, 'Spike Cut CutLeft must be greater or equal to Botm cut!');
    assert(P.botm.run==0 || P.botm.Tf <= P.spikeCutting.Tf, 'Spike Cut Tf must be greater or equal to Botm Tf!');
    assert(P.botm.run==0 ||P.spikeCutting.cutLeft - P.botm.cutLeft + P.botm.Tf <= P.spikeCutting.Tf, 'Spike Cut for Botm out of bounds!');
    assert(P.spikeCutting.cutLeft >= P.spikeAlignment.cutLeft, 'Cannot align further left than cut !');
    assert(P.spikeAlignment.cutLeft - P.featureExtraction.cutLeft >= 0, 'Connot prewhiten spikes with this cutleft!');
    % Save header
    readme = 'The files in this folder were created by the sort.m script.';
    
    % Define global names
    S.P = P;
    S.name = name;
    S.srate = DS.getSamplesPerSecond();
    S.nC = size(DS,2);
    S.dpath = dpath;
    S.dprefix = fullfile(dpath, name);
    
    if ~exist(S.dpath, 'file')
        fprintf('Output directory does not exists. Trying to create...')
        mkdir(S.dpath);
        fprintf(' success.\n')
    end    
    
    % Define file names
    S.files.artefact_file          = [S.dprefix '.010artefacts.mat'];
    S.files.spike_det_file         = [S.dprefix '.020spikes_det.mat'];
    S.files.spike_det_merged_file  = [S.dprefix '.030spikes_det_merged.mat'];
    S.files.spike_cut_file         = [S.dprefix '.040spikes_cut.mat'];
    S.files.spike_aligned_file     = [S.dprefix '.050spikes_aligned.mat'];
    S.files.cov_file               = [S.dprefix '.060cov.mat'];
    S.files.prewh_spike_file       = [S.dprefix '.070spikes_prewhitened.mat'];
    S.files.fet_spike_file         = [S.dprefix '.080spikes_features.mat'];
    S.files.meanshift_spike_file   = [S.dprefix '.090clusters_meanshift.mat'];
    S.files.botm_matching_file     = [S.dprefix '.100botm_matching.mat'];    
    S.files.merge_ms_clusters      = [S.dprefix '.110clusters_meanshift_merged.mat'];     
    S.files.botm_file              = [S.dprefix '.120botm.mat'];
    S.files.botm_templates_file    = [S.dprefix '.130botm_templates.mat'];
    S.files.resArtefact_file       = [S.dprefix '.140residual_artefacts.mat'];
    S.files.botm_templates_cleaned_file = [S.dprefix '.150botm_cleaned_templates.mat'];
    save([fullfile(dpath, name) '.P.mat'], 'S', 'readme');
    
    % Define constants
    nC = size(DS,2);    

    % RUN
    runArtefactDetection();
    runSpikeDetection();
    mergeDetectedSpikes();
    noiseEstimation();
    cutSpikes();
    alignSpikes();    
    prewhitenSpikes();
    fetExtraction();
    runMeanShift();
    runBOTMMatching();
    mergeClustersMS();
    runBOTM();
%     computeBOTMTemplates();
%     residualArtefactDetection();    
%     computeBOTMCleanedTemplates();
%     mergeClustersBotm();
   
    disp('############################');
    disp('All done.!');
    
    
    %----------------------------------------------------------------------
    % Detected sp:  +++++ ++ +++++++++++++ +++++++ ++++++++++ ++++++ +++++
    % Cut spikes :  + +++ ++ + + + ++ + ++ + + + + +++ ++ + + ++ +++ +++++
    % Aligned sp :  + ++     +   +      ++     + +      +   +  +  ++ ++  +
    % Prewhitend :  + ++     +   +      ++     + +      +   +  +  ++ ++  +
    % Features   :  + ++     +   +      ++     + +      +   +  +  ++ ++  +
    % Clustered  :     +     +           +     + +                 + +
    % Matched    :  + +++ ++ + + + ++ + ++ + + + + +++ ++ + + ++ +++ +++++
    % BOTM       :  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
    %% --------------------------------------------------------------------
    function runArtefactDetection()
        artefactDetection.epochs = [];        
        if P.artefactDetection.use
            if exist(S.files.artefact_file, 'file')
                disp('artefacts already detected');
                load(S.files.artefact_file, 'artefactDetection');
            else
                disp('detecting artefacts...');
                A = ana.moritzheimdahl.ArtefactDetector(DS, P.artefactDetection.width, P.artefactDetection.threshold);
        %         mysort.plot.SliderDataAxes({DS, A}, 'channelSpacers', [100 0]);
                artefactDetection.epochs = A.getEpochs();
                save(S.files.artefact_file, 'artefactDetection');
                disp('Done.')    
            end
        end
        S.artefactDetection = artefactDetection; clear artefactDetection;
    end

    %% --------------------------------------------------------------------
    function runSpikeDetection()
        if exist(S.files.spike_det_file, 'file')
            disp('spikes already detected');
            spikeDetection = [];
            load(S.files.spike_det_file);
        else
            disp('Detecting spikes...');
            L = P.spikeDetection.maxDataLength;
            spikeDetection.spikesDetectedUp = {};
            spikeDetection.pks_up = {};
            if ~isempty(strfind(P.spikeDetection.method, '+'))
                [spikeDetection.spikesDetectedUp spikeDetection.pks_up] = ...
                    DS.detectSpikes('Len', L, 'energyfun', @(x) x, 'minPeakDistance', P.spikeDetection.minDist, 'thr', P.spikeDetection.thr);
            end        
            spikeDetection.spikesDetectedDown = {};
            spikeDetection.pks_down = {};
            if ~isempty(strfind(P.spikeDetection.method, '-'))
                [spikeDetection.spikesDetectedDown spikeDetection.pks_down] = ...
                    DS.detectSpikes('Len', L, 'energyfun', @(x) -x, 'minPeakDistance', P.spikeDetection.minDist, 'thr', P.spikeDetection.thr);
            end
            save(S.files.spike_det_file, 'spikeDetection');
            disp('Done.')    
        end
        S.spikeDetection = spikeDetection; clear spikeDetection;
    end

    %% --------------------------------------------------------------------
    function mergeDetectedSpikes()
        if exist(S.files.spike_det_merged_file, 'file')
            disp('spikes already merged');
            spikeDetectionMerged = [];
            load(S.files.spike_det_merged_file);
        else
            disp('Merging spikes trains...');
            allspikes = [cell2mat(S.spikeDetection.spikesDetectedUp)   cell2mat(S.spikeDetection.pks_up);
                         cell2mat(S.spikeDetection.spikesDetectedDown) cell2mat(S.spikeDetection.pks_down)];
            nSp = size(allspikes,1);
            % add the channel information to allspikes
            allspikes(nSp,3) = 0;
            nextIdx = 1;
            LL = max(length(S.spikeDetection.spikesDetectedUp),...
                     length(S.spikeDetection.spikesDetectedDown));
            for i=1:LL
                nSpUp = 0;
                nSpDo = 0;
                if ~isempty(S.spikeDetection.spikesDetectedUp)
                    nSpUp = length(S.spikeDetection.spikesDetectedUp{i});
                end
                if ~isempty(S.spikeDetection.spikesDetectedDown)
                    nSpDo = length(S.spikeDetection.spikesDetectedDown{i});
                end
                allspikes(nextIdx:nextIdx+nSpUp+nSpDo-1,3) = i;
                nextIdx = nextIdx+nSpUp+nSpDo;
            end
            allspikes = sortrows(allspikes,1);
            fprintf('%d spikes found total before merging\n', nSp);
            % Check for oscillation groups and remove those before
            % artifact removal
            if P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude(1) > 0
                allspikes = mysort.spiketrain.findOscillationGroups(allspikes, ...
                    P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude(3), ...
                    P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude(1),...
                    P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude(2));
            end
            % remove spikes in artefact epochs
            removeIdx = mysort.epoch.findPointsInEpochs(allspikes(:,1), S.artefactDetection.epochs)>0;               
            allspikes(removeIdx,:) = [];
            allspikes(abs(allspikes(:,2))>P.spikeDetection.removeEventsWithAbsAmplitudeLargerThan,:) = [];
            fprintf('%d spikes after artefact removal (%d removed)\n', size(allspikes,1), nSp-size(allspikes,1));

            allspikes  = mysort.spiketrain.mergeSingleElectrodeDetectedSpikes(allspikes, P.spikeDetection.mergeSpikesMaxDist);
            % remove spikes at the borders
            allspikes(allspikes(:,1) < P.spikeCutting.Tf,:) = [];
            allspikes(allspikes(:,1) > size(DS,1),:) = [];
            spikeDetectionMerged.allspikes = allspikes;
            spikeDetectionMerged.ts = allspikes(:, 1);
            save(S.files.spike_det_merged_file, 'spikeDetectionMerged');
            disp('Done.')    
        end
        S.spikeDetectionMerged = spikeDetectionMerged; 
        fprintf('Found %d spikes after merging.\n', length(spikeDetectionMerged.ts));
    end

    %% --------------------------------------------------------------------
    function cutSpikes()
        if exist(S.files.spike_cut_file, 'file')
            disp('spikes already cut');
            spikeCut = [];
            load(S.files.spike_cut_file);
        else
            disp('Cutting spikes...');
            nSpikesDetected = length(S.spikeDetectionMerged.ts);
            if nSpikesDetected > P.spikeCutting.maxSpikes
                cutIdx = randperm(nSpikesDetected);
                spikeCut.cutIdx = sort(cutIdx(1:P.spikeCutting.maxSpikes));
            else
                spikeCut.cutIdx = 1:nSpikesDetected;
            end
            spikeCut.wfs = DS.getWaveform(S.spikeDetectionMerged.ts(spikeCut.cutIdx),...
                P.spikeCutting.cutLeft, P.spikeCutting.Tf);            
            save(S.files.spike_cut_file, 'spikeCut', '-v7.3');
            disp('Done.')    
        end
        S.spikeCut = spikeCut; clear spikeCut;
    end

    %% --------------------------------------------------------------------
    function noiseEstimation()
        if exist(S.files.cov_file, 'file')
            disp('cov already estimated...');
            noise = [];
            load(S.files.cov_file);
        else
            disp('Estimating Cov...');
            s1 = S.spikeDetectionMerged.ts-P.noiseEstimation.minDistFromSpikes;
            s2 = S.spikeDetectionMerged.ts+P.noiseEstimation.minDistFromSpikes;
            if isempty(P.spikeDetection.maxDataLength)
                L = size(DS,1);
            else
                L  = P.spikeDetection.maxDataLength;
            end
            noise.epochs = mysort.epoch.flip(mysort.epoch.merge([s1 s2]), L);
            LN = P.noiseEstimation.minLength;
            maxTf = max([P.botm.Tf P.featureExtraction.Tf P.spikeAlignment.Tf]); %P.spikeCutting.Tf
            Cest = mysort.noise.Covest2(DS, 'maxLag', maxTf,...
                'maxSamples', LN, 'noiseEpochs', noise.epochs, 'forceMethod', 'xcorr');
            noise.C_time_cut = mysort.noise.ccol2Cte(Cest.CCol, maxTf);
            noise.C_time_aligned = mysort.noise.ccol2Cte(Cest.CCol, P.spikeAlignment.Tf);
            noise.meanNoiseStd = sqrt(mean(diag(noise.C_time_cut)));
            noise.CestS = Cest.toStruct();    
            save(S.files.cov_file, 'noise');
            disp('Done.')    
        end
        S.noise = noise; clear noise;
    end

    %% --------------------------------------------------------------------
    function alignSpikes()
        if ~isempty(P.spikeAlignment.method); 
            if exist(S.files.spike_aligned_file, 'file')
                disp('spikes already aligned');
                spikeAligned = [];
                load(S.files.spike_aligned_file);
            else
                disp('Aligning spikes...');
                idxstart = 1 + P.spikeCutting.cutLeft - P.spikeAlignment.cutLeft;
                idxstop  = idxstart + P.spikeAlignment.Tf - 1;
                
                nSpCut = size(S.spikeCut.wfs,1);
                if P.spikeAlignment.maxSpikes < nSpCut
                    alignIdx = randperm(nSpCut);
                    spikeAligned.alignIdx = sort(alignIdx(1:P.spikeAlignment.maxSpikes));
                else
                    spikeAligned.alignIdx = 1:nSpCut;
                end
                spikeAligned.unalignedwfs = mysort.wf.vSubsel(S.spikeCut.wfs(spikeAligned.alignIdx,:),...
                    nC, idxstart:idxstop);    

                spikeAligned.tau = zeros(size(spikeAligned.unalignedwfs,1),1);
                spikeAligned.maxIdx = P.spikeAlignment.cutLeft;
                
                spikeAligned.restrictToIdx = spikeAligned.maxIdx-2:spikeAligned.maxIdx+5;
                tic
                if strcmp(P.spikeAlignment.method, 'onMax')
                    [spikeAligned.wfs spikeAligned.tau] = ...
                        mysort.wf.vAlignOnMax(spikeAligned.unalignedwfs, nC,...
                        'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                elseif strcmp(P.spikeAlignment.method, 'onUpsampledMean')
                    [spikeAligned.wfs spikeAligned.tau] = ...
                        mysort.wf.vAlignOnUpsampleMean(spikeAligned.unalignedwfs, nC,...
                        'maxIdx', P.spikeAlignment.maxIdx,...
                        'maxIter', P.spikeAlignment.maxIterations, ...
                        'initAlignment', P.spikeAlignment.initAlignment);   
                elseif strcmp(P.spikeAlignment.method, 'onUpsampledMax')
                    [spikeAligned.tau spikeAligned.wfs] = ...
                        mysort.util.alignWaveformsUpsampleMax(spikeAligned.unalignedwfs, nC,...
                        'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx, 'nIter', 2);
                elseif strcmp(P.spikeAlignment.method, 'onUpsampledMin')
                    [spikeAligned.tau spikeAligned.wfs] = ...
                        mysort.wf.alignWaveformsUpsampleMax(-spikeAligned.unalignedwfs, nC,...
                        'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx, 'nIter', 2); 
                    spikeAligned.wfs = -spikeAligned.unalignedwfs;
                elseif strcmp(P.spikeAlignment.method, 'onMin')
                    [spikeAligned.wfs spikeAligned.tau] = ...
                        mysort.wf.vAlignOnMax(-spikeAligned.unalignedwfs, nC,...
                        'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                    spikeAligned.wfs = -spikeAligned.wfs;
                elseif strcmp(P.spikeAlignment.method, 'onAverageMax')
                    [spikeAligned.tau spikeAligned.wfs] = ...
                        mysort.wf.vAlignOnAverageMaxSample(spikeAligned.unalignedwfs, nC,...
                        'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                elseif strcmp(P.spikeAlignment.method, 'onAverageMin')
                    [spikeAligned.tau spikeAligned.wfs] = ...
                        mysort.wf.vAlignOnAverageMaxSample(-spikeAligned.unalignedwfs, nC,...
                        'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                    spikeAligned.wfs = -spikeAligned.wfs;
                elseif strcmp(P.spikeAlignment.method, 'none')
                    spikeAligned.wfs = spikeAligned.unalignedwfs;
                else
                    error(['unkown alignement method' P.spikeAlignment.method]);
                end
                toc
                % cut away alignement artefacts at the ends
%                 spikeAligned.unalignedwfs([1:max(1,max(spikeAligned.tau)) end+min(0,min(spikeAligned.tau))+1:end],:) = [];
                save(S.files.spike_aligned_file, 'spikeAligned', '-v7.3');
                disp('Done.')
            end
        end
        S.spikeAligned = spikeAligned; clear spikeAligned;        
    end

    %% --------------------------------------------------------------------
    function prewhitenSpikes()
        if exist(S.files.prewh_spike_file, 'file')
            disp('spikes already prewhitened...');
            spikePrewhitened = [];
            load(S.files.prewh_spike_file);
        else
            disp('prewhitening spikes...');
            % Cut the aligned waveforms to get rid of alignement artefacts
            % and reduce to the final waveform that will be used for
            % feature extraction
            idxstart = 1 + P.spikeAlignment.cutLeft - P.featureExtraction.cutLeft;
            idxstop  = idxstart + P.featureExtraction.Tf - 1;
            spikePrewhitened.wfs = mysort.wf.vSubsel(S.spikeAligned.wfs, nC, idxstart:idxstop);      
            % Build the noise covariance matrix and load it
            ccol_loaded = S.noise.CestS.CCol;
            ccol_loaded(1:nC, 1:nC) = ccol_loaded(1:nC, 1:nC) + diag(diag(ccol_loaded(1:nC, 1:nC)));
            spikePrewhitened.ccol_loaded = ccol_loaded/2;
            spikePrewhitened.C = mysort.noise.ccol2Cte(spikePrewhitened.ccol_loaded, P.featureExtraction.Tf);
            spikePrewhitened.U = chol(spikePrewhitened.C);
            % Prewhiten
            spikePrewhitened.wfs = spikePrewhitened.wfs/spikePrewhitened.U; 
            save(S.files.prewh_spike_file, 'spikePrewhitened', '-v7.3');
            disp('Done.')    
        end
        S.spikePrewhitened = spikePrewhitened; clear spikePrewhitened;
    end

    %% --------------------------------------------------------------------
    function fetExtraction()
        if exist(S.files.fet_spike_file , 'file')
            disp('features already calculated...');
            spikeFeatures = [];
            load(S.files.fet_spike_file );
        else
            disp('Calculating features...');
            spikeFeatures.X = mysort.util.dimReductionPCA(S.spikePrewhitened.wfs,...
                P.featureExtraction.nDims, [], 3*1000000);
            save(S.files.fet_spike_file, 'spikeFeatures', '-v7.3');
            disp('Done.')    
        end
        S.spikeFeatures = spikeFeatures; clear spikeFeatures;
    end
    
    %% --------------------------------------------------------------------
    function runMeanShift()
        if exist(S.files.meanshift_spike_file, 'file')
            disp('Already clustered...');
            clustering = [];
            load(S.files.meanshift_spike_file);
        else
            disp('Mean Shift Clustering...');
            if ~isempty(P.clustering.meanShiftBandWidth)
                clustering.bandwidth = P.clustering.meanShiftBandWidth;
            else
                clustering.bandwidth = sqrt(1.1*P.featureExtraction.nDims);
            end
            nSpFet = size(S.spikeFeatures.X,1);
            if P.clustering.maxSpikes < nSpFet
                clusterIdx = randperm(nSpFet);
                clustering.clusterIdx = sort(clusterIdx(1:P.clustering.maxSpikes));
            else
                clustering.clusterIdx = 1:nSpFet;
            end
            X = S.spikeFeatures.X(clustering.clusterIdx,:)';
            S.maxNBandwidthIncreases = 1;
            S.bandwidthIncreaseFactor = 1.3;
            tic
%             [clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(X, clustering.bandwidth);
            [clustCent,point2cluster,clustMembsCell] = MeanShiftClusterIncreaseBW(X, clustering.bandwidth, 0 , P.clustering.minSpikesPerCluster, S.maxNBandwidthIncreases, S.bandwidthIncreaseFactor)
            toc  
            [clustering.ids clustering.clusterCenter] = MeanShiftClusterBundleResult(X', clustMembsCell, P.clustering.minSpikesPerCluster);
            clustering.classes = unique(clustering.ids);
            clusteredAlignedSpikes = S.spikeAligned.wfs(clustering.clusterIdx,:);
            clusteredCutSpikes = S.spikeCut.wfs(S.spikeAligned.alignIdx(clustering.clusterIdx),:);
            
            clustering.templatesCut = mysort.util.calculateClassMeans(clusteredCutSpikes, clustering.ids, 'usemedian');
            clustering.templatesAligned = mysort.util.calculateClassMeans(clusteredAlignedSpikes, clustering.ids, 'usemedian');
            save(S.files.meanshift_spike_file, 'clustering', '-v7.3');
        end
        S.clustering = clustering; clear clustering;
    end
    %% --------------------------------------------------------------------
    function runBOTMMatching()
        if exist(S.files.botm_matching_file, 'file')
            disp('botm matched with lda...');
            clusteringMatched = [];
            load(S.files.botm_matching_file);
        else
            disp('BOTM matching...');
            nSpCut = size(S.spikeCut.wfs,1);
            clusteringMatched.ids = zeros(nSpCut,1);
            clusteringMatched.template_ids = 0;
            clusteringMatched.spikeCutAligned = [];
            clusteringMatched.templates = [];
            clusteringMatched.maxTausPerSpike = zeros(nSpCut,1);
            clusteringMatched.ts = S.spikeDetectionMerged.ts(S.spikeCut.cutIdx);
                        
            if size(S.clustering.templatesCut,1) < 2
                disp('Cannot run BOTM Matching, not enough templates!');
            else
                T = S.clustering.templatesAligned; % leave a noise template
                nT = size(T,1);
                C = (S.noise.C_time_aligned +diag(diag(S.noise.C_time_aligned)))/2;

                % Do matching only on temporal window defined by the aligned
                % templates
                idxstart = 1 + P.spikeCutting.cutLeft - P.spikeAlignment.cutLeft;
%                 idxstop  = idxstart + P.spikeAlignment.Tf - 1;

                [D maxTaus TupShiftDown FupShiftDown tauRange tauIdx subTauRange subTauIdx] = mysort.wf.vTemplateMatching(...
                    S.spikeCut.wfs, T, nC, idxstart, 'maxShift', 5, 'upsample', 5, 'noiseCovariance', C);
                F = FupShiftDown(:,:,1);
                EF = diag(T*F');                     % compute energies
                Prior  = P.templateMatchingCut.prior;                
                DISCR = D - .5 * repmat(EF', nSpCut, 1)  + log(Prior);  % compute botm dicriminant
                clusteringMatched.maxTausPerSpikeAndFilter = maxTaus;
                clusteringMatched.tauRange = tauRange;
                clusteringMatched.subTauRange = subTauRange;
                clusteringMatched.tauIdx = tauIdx;
                clusteringMatched.subTauIdx = subTauIdx;
                clusteringMatched.TupShiftDown = TupShiftDown;
                clusteringMatched.FupShiftDown = FupShiftDown;
                
                [MaxD Didx] = max(DISCR,[],2);
                uDidx = unique(Didx);
                for k = 1:length(uDidx)
                    locDidx = Didx == uDidx(k);
                    clusteringMatched.maxTausPerSpike(locDidx,1) = maxTaus(locDidx,uDidx(k));
                    clusteringMatched.maxTausIdxPerSpike(locDidx,1) = tauIdx(locDidx,uDidx(k));
                    clusteringMatched.maxSubTausIdxPerSpike(locDidx,1) = subTauIdx(locDidx,uDidx(k));
                end
                notNoiseIdx = MaxD>0;
                clusteringMatched.ids(notNoiseIdx,1) = Didx(notNoiseIdx);
                clusteringMatched.template_ids = unique(clusteringMatched.ids);

                if 0 %~isempty(P.templateMatchingCut.residualPeakSmallerThanStdNoise)
                    % compute the "explained energy"
                    notNoiseIdx = find(notNoiseIdx);
                    nSmatched = length(notNoiseIdx);
                    E = zeros(size(S.spikeCut.wfs));
                    tTf = size(TupShiftDown,2);
                    for s = 1:8 %nSmatched
                        myId = clusteringMatched.ids(notNoiseIdx(s));
                        mySubTauIdx = clusteringMatched.maxSubTausIdxPerSpike(notNoiseIdx(s));
                        myTau = floor(clusteringMatched.maxTausPerSpike(notNoiseIdx(s)));
                        myStartIdx = idxstart-myTau;
                        E(notNoiseIdx(s),myStartIdx:myStartIdx+tTf-1) = squeeze(TupShiftDown(myId,:,mySubTauIdx));
                    end
                    R = S.spikeCut.wfs-E;
                    figure;
                    subplot(2,1,1);
                    plot(S.spikeCut.wfs(1:4,:)');
                    hold on
                    plot(E(1:4,:)');
                    subplot(2,1,2); plot(R(1:4,:)');
                end
                
                % shift the spikes to compute the templates
                clusteringMatched.spikeCutAligned = mysort.wf.vShift(S.spikeCut.wfs, nC, -round(clusteringMatched.maxTausPerSpike), 1);
                clusteringMatched.templates = mysort.util.calculateClassMeans(...
                    clusteringMatched.spikeCutAligned, clusteringMatched.ids, 'usemedian');  
                clusteringMatched.templatesIDs = unique(clusteringMatched.ids);
                clusteringMatched.ts = S.spikeDetectionMerged.ts(S.spikeCut.cutIdx)+clusteringMatched.maxTausPerSpike;
                if 0
                    N = 1000;
                    figure; subplot(2,1,1); plot(S.spikeCut.wfs(1:N,:)');
                    subplot(2,1,2);
                    plot(clusteringMatched.spikeCutAligned(1:N,:)');
                end
            end
            save(S.files.botm_matching_file, 'clusteringMatched', '-v7.3'); 
            disp('Done.')             
        end
        S.clusteringMatched = clusteringMatched;
    end
    %% --------------------------------------------------------------------
    function mergeClustersMS()
        if ~P.mergeTemplates.merge
            return
        end        
        if exist(S.files.merge_ms_clusters, 'file')
            disp('Already merged...');
            clusteringMerged = [];
            load(S.files.merge_ms_clusters);
        else   
            disp('Merging MS clusters...');
            
            % in S.clusteringMatched.ids the ids are going from 0 to nT with
            % 0 being not matched and noise, thus there is one "real" template less
            % than units (if any spike has an ID of 0 !!)
            realTemplateIdx = find(S.clusteringMatched.template_ids>0);
            nT = length(realTemplateIdx);
            if nT == 0
                % we have only the noise cluster
                D = [];
                maxT = [];
                groups = {};
            else
                % resample only templates that do not have id==0
                resampledTemplates = mysort.util.resampleTensor(mysort.wf.v2t(...
                    S.clusteringMatched.templates(realTemplateIdx,:), nC), 3, 1);
%                 mysort.plot.waveforms(resampledTemplates, 'IDs', S.clusteringMatched.template_ids(realTemplateIdx), 'stacked', 0);

                % the indices in the groups are indices into
                % realTemplateIdx which is index into S.clusteringMatched.template_ids
                [groups maxT D] = ana.mergeTemplates(resampledTemplates, S.noise.meanNoiseStd,...
                    'maxRelativeDistance', P.mergeTemplates.ifMaxRelDistSmallerPercent/100,...
                    'minCorrelation', P.mergeTemplates.atCorrelation);
            end
            clusteringMerged.D = D;
            clusteringMerged.maxT = maxT;
            clusteringMerged.groups = groups;            
            clusteringMerged.templates = zeros(length(groups), P.spikeCutting.Tf*nC);
            clusteringMerged.ids = zeros(length(S.clusteringMatched.ids),1);

            % we need to merge now the indices in ids accoring to the groups in
            % groups. But groups points into realTemplateIdx !!
            for i=1:length(groups)
                group_ids = S.clusteringMatched.template_ids(realTemplateIdx(clusteringMerged.groups{i}));
                idx = ismember(S.clusteringMatched.ids, group_ids);
                clusteringMerged.ids(idx) = i;
                clusteringMerged.templates(i,:) = median(S.clusteringMatched.spikeCutAligned(idx,:),1);
            end 
%             mysort.plot.waveforms(S.spikeAligned.wfs, 'nC', nC, 'IDs', clusteringMerged.ids, 'stacked', 0)
            
            save(S.files.merge_ms_clusters, 'clusteringMerged', '-v7.3');
            fprintf('Templates after merging: %d\n', length(clusteringMerged.groups));
        end
        S.clusteringMerged = clusteringMerged;
    end    

    %% --------------------------------------------------------------------
    function runBOTM()
        if ~P.botm.run
            return
        end

        if exist(S.files.botm_file, 'file')
            disp('already sorted with botm...');
            botm = [];
            load(S.files.botm_file);
        else
            disp('BOTM sorting...');
            T = S.clusteringMerged.templates;
            Cest = S.noise.CestS;
            Cest.CCol(1:nC,:) = Cest.CCol(1:nC,:) + diag(diag(Cest.CCol(1:nC,:)));
            Cest.CCol = Cest.CCol/2;
            Cest = mysort.noise.Covest2(Cest);
            idxstart = 1 + P.spikeCutting.cutLeft - P.botm.cutLeft;
            idxstop  = idxstart + P.botm.Tf - 1;
            botm.templates = mysort.wf.vSubsel(T, nC, idxstart:idxstop);
            botmsorter = mysort.sorters.BOTM(Cest, P.botm.Tf, botm.templates,...
                        'upsample', 3, 'spikePrior', P.botm.prior,...
                        'max_num_of_channels_per_template', 30, 'adaptOnInit', 1);    
            
            botmsorter.DH = DS;    
            botm.gdfUnclean = botmsorter.sort(DS);
            % Remove artefacts from gdf
            botm.removeIdx = mysort.epoch.findPointsInEpochs(botm.gdfUnclean(:,2), S.artefactDetection.epochs)>0;               
            botm.gdf = botm.gdfUnclean(~botm.removeIdx,:);
            save(S.files.botm_file, 'botm');
            disp('Done.')    
        end
        S.botm = botm; 
    end
    %% --------------------------------------------------------------------
    function computeBOTMTemplates()
        if ~isempty(gdfbotm)
            if exist(S.files.botm_templates_file, 'file')
                disp('loading botm templates...');
                D = load(S.files.botm_templates_file);
                templatesAfterBOTM = D.templatesAfterBOTM;
                clear D
            else
                disp('computing botm templates...');
                uclasses = unique(gdfbotm(:,1));
                templatesAfterBOTM = zeros(length(uclasses), P.spike_cut_Tf*nC);
                for i=1:length(uclasses)
                    myts = gdfbotm(gdfbotm(:,1)==uclasses(i),2);
                    if size(myts,1) > 1000
                        myts = myts(1:1000);
                    end
                    spikes = DS.getWaveform(myts, P.spike_cut_cutLeft, P.spike_cut_Tf);
                    templatesAfterBOTM(i,:) = mean(spikes,1);
                end
                save(S.files.botm_templates_file, 'templatesAfterBOTM');
                disp('Done.')    
            end            
        end
        S.botmDetails = botmDetails; clear botmDetails;
    end
    %% --------------------------------------------------------------------
    function residualArtefactDetection()
        if ~isempty(gdfbotm)
            if exist(S.files.resArtefact_file, 'file')
                disp('residual artefacts already detected');
                D = load(S.files.resArtefact_file);
                gdfbotmcleaned = D.gdfbotmcleaned;
                resArtefactEpochs = D.resArtefactEpochs; clear D;
            else
                disp('Residual atrefact detection ...');
                SpSoBOTM = mysort.spiketrain.SpikeSortingContainer('botm', gdfbotm, ...
                    'templateCutLeft', S.P.spike_cut_cutLeft, ...
                    'templateCutLength', S.P.spike_cut_Tf, ...
                    'templateWfs', mysort.wf.v2t(templatesAfterBOTM, nC),...
                    'wfDataSource', DS, 'nMaxSpikesForTemplateCalc', 1000);
                DS.addSpikeSorting(SpSoBOTM); DS.setActiveSpikeSortingIdx('botm')            
                DS.bReturnSortingResiduals = 1;
                A = ana.moritzheimdahl.ResidualArtefactDetector(DS);
%                 mysort.plot.SliderDataAxes({DS, A}, 'channelSpacers', [5000 0]);
                resArtefactEpochs = A.getEpochs();
                gdfbotmcleaned = gdfbotm;
                removeIdx = mysort.epoch.findPointsInEpochs(gdfbotm(:,2), resArtefactEpochs)>0;
                gdfbotmcleaned(removeIdx,:) = [];      
                save(S.files.resArtefact_file, 'resArtefactEpochs', 'gdfbotmcleaned');
                DS.bReturnSortingResiduals = 0;
                disp('Done.')    
            end
        end
        S.botmResiduals = botmResiduals; clear botmResiduals;
    end
    %% --------------------------------------------------------------------
    function computeBOTMCleanedTemplates()
        if ~isempty(gdfbotmcleaned)
            if exist(S.files.botm_templates_cleaned_file, 'file')
                disp('loading botm cleaned templates...');
                D = load(S.files.botm_templates_cleaned_file);
                templatesAfterBOTMcleaned = D.templatesAfterBOTMcleaned;
                clear D
            else
                disp('computing botm cleaned templates...');
                uclasses = unique(gdfbotmcleaned(:,1));
                templatesAfterBOTMcleaned = zeros(length(uclasses), P.spike_cut_Tf*nC);
                for i=1:length(uclasses)
                    myts = gdfbotmcleaned(gdfbotmcleaned(:,1)==uclasses(i),2);
                    if size(myts,1) > 1000
                        myts = myts(1:1000);
                    end
                    spikes = DS.getWaveform(myts, P.spike_cut_cutLeft, P.spike_cut_Tf);
                    templatesAfterBOTMcleaned(i,:) = mean(spikes,1);
                end
                save(S.files.botm_templates_cleaned_file, 'templatesAfterBOTMcleaned');
                disp('Done.')    
            end            
        end
        S.botmCleaned = botmCleaned; clear botmCleaned;                
    end
%     %% --------------------------------------------------------------------
%     function mergeClustersBotm()
%         if 1 %isempty(gdfbotmcleaned)
%             return
%         end
%     end
end