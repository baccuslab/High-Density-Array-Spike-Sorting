classdef sortDS < handle
    properties (SetAccess=private)
        
    end
    properties
        DS
        P
        files
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = sortDS(DS, dpath, varargin)  
            % Data Specific
            self.P.channelIdx = 1:size(DS,2);
            
            % Artefact Detection
            self.P.artefactDetection.use = 0;
            self.P.artefactDetection.width = 350;
            self.P.artefactDetection.threshold = size(DS,2)*41; % should be sampling frequency dependent!

            % Spike Detection
            self.P.spikeDetection.method = '+-';
            self.P.spikeDetection.thr = 3.5;
            self.P.spikeDetection.minDist = 25;
            self.P.spikeDetection.maxDataLength = [];
            self.P.spikeDetection.mergeSpikesMaxDist = 18; % do not go below 13 for HDMEA!
            self.P.spikeDetection.removeEventsWithAbsAmplitudeLargerThan = 15000;
            self.P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude = [0 0 0]; % should remove atrifacts that are oscillations (successive detections of peaks with alternating signs). First parameter is N, second the minHeight, third the maxDist within a group

            % Spike cutting
            self.P.spikeCutting.maxSpikes = 100000;
            self.P.spikeCutting.Tf = 55;
            self.P.spikeCutting.cutLeft = 10;

            % Spike Alignment
            self.P.spikeAlignment.method = 'onUpsampledMean';
            self.P.spikeAlignment.maxSpikes = 60000;
            self.P.spikeAlignment.Tf = 25;
            self.P.spikeAlignment.cutLeft = 6;
            self.P.spikeAlignment.initAlignment = '-';
            self.P.spikeAlignment.maxIdx = self.P.spikeAlignment.cutLeft + 1;

            self.P.spikeAlignment.maxIterations = 30;

            % Noise estimation    
            self.P.noiseEstimation.minLength = 2000000;
            self.P.noiseEstimation.minDistFromSpikes = 60;

            % Feature extraction
            self.P.featureExtraction.Tf = 15;
            self.P.featureExtraction.cutLeft = 3;
            self.P.featureExtraction.nDims = 6;

            % Clustering
            self.P.clustering.maxSpikes = 30000;
            self.P.clustering.meanShiftBandWidth = []; % default is sqrt(1.1*nDims)
        %     self.P.clustering.bandwidthMultMaxDistSpeedMS = 5;
            self.P.clustering.minSpikesPerCluster = 10;

            % Template Matching on Cut Spikes
            self.P.templateMatchingCut.prior = .0001;
            self.P.templateMatchingCut.residualPeakSmallerThanStdNoise = [];

            % Merge Templates after template matching
            self.P.mergeTemplates.merge = 1;
            self.P.mergeTemplates.upsampleFactor = 3;
            self.P.mergeTemplates.atCorrelation = .96;
        %     self.P.mergeTemplates.ifMeanDistSmaller = .1;    % in units std of noise
            self.P.mergeTemplates.ifMaxDistSmaller = 2.5;      % in units std of noise
            self.P.mergeTemplates.ifMaxRelDistSmallerPercent = 25;

            % Run BOTM Template Matching on whole data
            self.P.botm.run = 0;
            self.P.botm.Tf = 55;
            self.P.botm.cutLeft = 10;
            self.P.botm.prior = .0001;

            self.P = mysort.util.parseInputs(self.P, varargin, 'error');

            % Checks
            assert(self.P.botm.cutLeft <= self.P.spikeCutting.cutLeft, 'Spike Cut CutLeft must be greater or equal to Botm cut!');
            assert(self.P.botm.Tf      <= self.P.spikeCutting.Tf, 'Spike Cut Tf must be greater or equal to Botm Tf!');
            assert(self.P.spikeCutting.cutLeft - self.P.botm.cutLeft + self.P.botm.Tf <= self.P.spikeCutting.Tf, 'Spike Cut for Botm out of bounds!');
            assert(self.P.spikeCutting.cutLeft >= self.P.spikeAlignment.cutLeft, 'Cannot align further left than cut !');
            assert(self.P.spikeAlignment.cutLeft - self.P.featureExtraction.cutLeft >= 0, 'Connot prewhiten spikes with this cutleft!');
            % Save header
            readme = 'The files in this folder were created by the sort.m script.';

            % Define global names
            self.name = name;
            self.srate = self.DS.getSamplesPerSecond();
            self.nC = size(DS,2);
            self.dpath = dpath;
            self.dprefix = fullfile(dpath, name);

            if ~exist(self.dpath, 'file')
                fprintf('Output directory does not exists. Trying to create...')
                mkdir(self.dpath);
                fprintf(' success.\n')
            end    

            % Define file names
            self.files.artefact_file          = [self.dprefix '.010artefacts.mat'];
            self.files.spike_det_file         = [self.dprefix '.020spikes_det.mat'];
            self.files.spike_det_merged_file  = [self.dprefix '.030spikes_det_merged.mat'];
            self.files.spike_cut_file         = [self.dprefix '.040spikes_cut.mat'];
            self.files.spike_aligned_file     = [self.dprefix '.050spikes_aligned.mat'];
            self.files.cov_file               = [self.dprefix '.060cov.mat'];
            self.files.prewh_spike_file       = [self.dprefix '.070spikes_prewhitened.mat'];
            self.files.fet_spike_file         = [self.dprefix '.080spikes_features.mat'];
            self.files.meanshift_spike_file   = [self.dprefix '.090clusters_meanshift.mat'];
            self.files.botm_matching_file     = [self.dprefix '.100botm_matching.mat'];    
            self.files.merge_ms_clusters      = [self.dprefix '.110clusters_meanshift_merged.mat'];     
            self.files.botm_file              = [self.dprefix '.120botm.mat'];
            self.files.botm_templates_file    = [self.dprefix '.130botm_templates.mat'];
            self.files.resArtefact_file       = [self.dprefix '.140residual_artefacts.mat'];
            self.files.botm_templates_cleaned_file = [self.dprefix '.150botm_cleaned_templates.mat'];
            save([fullfile(dpath, name) '.self.P.mat'], 'S', 'readme');

            % Define constants
            nC = size(DS,2);    

            % RUN
            self.runArtefactDetection();
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

        end
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
        function runArtefactDetection(self)
            artefactDetection.epochs = [];        
            if self.P.artefactDetection.use
                if exist(self.files.artefact_file, 'file')
                    disp('artefacts already detected');
                    load(self.files.artefact_file, 'artefactDetection');
                else
                    disp('detecting artefacts...');
                    A = ana.moritzheimdahl.ArtefactDetector(DS, self.P.artefactDetection.width, self.P.artefactDetection.threshold);
            %         mysort.plot.SliderDataAxes({DS, A}, 'channelSpacers', [100 0]);
                    artefactDetection.epochs = A.getEpochs();
                    save(self.files.artefact_file, 'artefactDetection');
                    disp('Done.')    
                end
            end
            self.artefactDetection = artefactDetection; clear artefactDetection;
        end

        %% --------------------------------------------------------------------
        function runSpikeDetection(self)
            if exist(self.files.spike_det_file, 'file')
                disp('spikes already detected');
                spikeDetection = [];
                load(self.files.spike_det_file);
            else
                disp('Detecting spikes...');
                L = self.P.spikeDetection.maxDataLength;
                spikeDetection.spikesDetectedUp = {};
                spikeDetection.pks_up = {};
                if ~isempty(strfind(self.P.spikeDetection.method, '+'))
                    [spikeDetection.spikesDetectedUp spikeDetection.pks_up] = ...
                        self.DS.detectSpikes('Len', L, 'energyfun', @(x) x, 'minPeakDistance', self.P.spikeDetection.minDist, 'thr', self.P.spikeDetection.thr);
                end        
                spikeDetection.spikesDetectedDown = {};
                spikeDetection.pks_down = {};
                if ~isempty(strfind(self.P.spikeDetection.method, '-'))
                    [spikeDetection.spikesDetectedDown spikeDetection.pks_down] = ...
                        self.DS.detectSpikes('Len', L, 'energyfun', @(x) -x, 'minPeakDistance', self.P.spikeDetection.minDist, 'thr', self.P.spikeDetection.thr);
                end
                save(self.files.spike_det_file, 'spikeDetection');
                disp('Done.')    
            end
            self.spikeDetection = spikeDetection; clear spikeDetection;
        end

        %% --------------------------------------------------------------------
        function mergeDetectedSpikes(self)
            if exist(self.files.spike_det_merged_file, 'file')
                disp('spikes already merged');
                spikeDetectionMerged = [];
                load(self.files.spike_det_merged_file);
            else
                disp('Merging spikes trains...');
                allspikes = [cell2mat(self.spikeDetection.spikesDetectedUp)   cell2mat(self.spikeDetection.pks_up);
                             cell2mat(self.spikeDetection.spikesDetectedDown) cell2mat(self.spikeDetection.pks_down)];
                nSp = size(allspikes,1);
                % add the channel information to allspikes
                allspikes(nSp,3) = 0;
                nextIdx = 1;
                for i=1:length(self.spikeDetection.spikesDetectedUp)
                    nSpUp = length(self.spikeDetection.spikesDetectedUp{i});
                    nSpDo = length(self.spikeDetection.spikesDetectedDown{i});
                    allspikes(nextIdx:nextIdx+nSpUp+nSpDo-1,3) = i;
                    nextIdx = nextIdx+nSpUp+nSpDo;
                end
                allspikes = sortrows(allspikes,1);
                fprintf('%d spikes found total before merging\n', nSp);
                % Check for oscillation groups and remove those before
                % artifact removal
                if self.P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude(1) > 0
                    allspikes = mysort.spiketrain.findOscillationGroups(allspikes, ...
                        self.P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude(3), ...
                        self.P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude(1),...
                        self.P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude(2));
                end
                % remove spikes in artefact epochs
                removeIdx = mysort.epoch.findPointsInEpochs(allspikes(:,1), self.artefactDetection.epochs)>0;               
                allspikes(removeIdx,:) = [];
                allspikes(abs(allspikes(:,2))>self.P.spikeDetection.removeEventsWithAbsAmplitudeLargerThan,:) = [];
                fprintf('%d spikes after artefact removal (%d removed)\n', size(allspikes,1), nSp-size(allspikes,1));

                allspikes  = mysort.spiketrain.mergeSingleElectrodeDetectedSpikes(allspikes, self.P.spikeDetection.mergeSpikesMaxDist);
                spikeDetectionMerged.allspikes = allspikes;
                spikeDetectionMerged.ts = allspikes(:, 1);
                save(self.files.spike_det_merged_file, 'spikeDetectionMerged');
                disp('Done.')    
            end
            self.spikeDetectionMerged = spikeDetectionMerged; 
            fprintf('Found %d spikes after merging.\n', length(spikeDetectionMerged.ts));
        end

        %% --------------------------------------------------------------------
        function cutSpikes(self)
            if exist(self.files.spike_cut_file, 'file')
                disp('spikes already cut');
                spikeCut = [];
                load(self.files.spike_cut_file);
            else
                disp('Cutting spikes...');
                nSpikesDetected = length(self.spikeDetectionMerged.ts);
                if nSpikesDetected > self.P.spikeCutting.maxSpikes
                    cutIdx = randperm(nSpikesDetected);
                    spikeCut.cutIdx = sort(cutIdx(1:self.P.spikeCutting.maxSpikes));
                else
                    spikeCut.cutIdx = 1:nSpikesDetected;
                end
                spikeCut.wfs = self.DS.getWaveform(self.spikeDetectionMerged.ts(spikeCut.cutIdx),...
                    self.P.spikeCutting.cutLeft, self.P.spikeCutting.Tf);            
                save(self.files.spike_cut_file, 'spikeCut', '-v7.3');
                disp('Done.')    
            end
            self.spikeCut = spikeCut; clear spikeCut;
        end

        %% --------------------------------------------------------------------
        function noiseEstimation(self)
            if exist(self.files.cov_file, 'file')
                disp('cov already estimated...');
                noise = [];
                load(self.files.cov_file);
            else
                disp('Estimating Cov...');
                s1 = self.spikeDetectionMerged.ts-self.P.noiseEstimation.minDistFromSpikes;
                s2 = self.spikeDetectionMerged.ts+self.P.noiseEstimation.minDistFromSpikes;
                if isempty(self.P.spikeDetection.maxDataLength)
                    L = size(DS,1);
                else
                    L  = self.P.spikeDetection.maxDataLength;
                end
                noise.epochs = mysort.epoch.flip(mysort.epoch.merge([s1 s2]), L);
                LN = self.P.noiseEstimation.minLength;
                maxTf = max([self.P.botm.Tf self.P.featureExtraction.Tf self.P.spikeAlignment.Tf]); %self.P.spikeCutting.Tf
                Cest = mysort.noise.Covest2(DS, 'maxLag', maxTf,...
                    'maxSamples', LN, 'noiseEpochs', noise.epochs, 'forceMethod', 'xcorr');
                noise.C_time_cut = mysort.noise.ccol2Cte(Cest.CCol, maxTf);
                noise.C_time_aligned = mysort.noise.ccol2Cte(Cest.CCol, self.P.spikeAlignment.Tf);
                noise.meanNoiseStd = sqrt(mean(diag(noise.C_time_cut)));
                noise.CestS = Cest.toStruct();    
                save(self.files.cov_file, 'noise');
                disp('Done.')    
            end
            self.noise = noise; clear noise;
        end

        %% --------------------------------------------------------------------
        function alignSpikes(self)
            if ~isempty(self.P.spikeAlignment.method); 
                if exist(self.files.spike_aligned_file, 'file')
                    disp('spikes already aligned');
                    spikeAligned = [];
                    load(self.files.spike_aligned_file);
                else
                    disp('Aligning spikes...');
                    idxstart = 1 + self.P.spikeCutting.cutLeft - self.P.spikeAlignment.cutLeft;
                    idxstop  = idxstart + self.P.spikeAlignment.Tf - 1;

                    nSpCut = size(self.spikeCut.wfs,1);
                    if self.P.spikeAlignment.maxSpikes < nSpCut
                        alignIdx = randperm(nSpCut);
                        spikeAligned.alignIdx = sort(alignIdx(1:self.P.spikeAlignment.maxSpikes));
                    else
                        spikeAligned.alignIdx = 1:nSpCut;
                    end
                    spikeAligned.unalignedwfs = mysort.wf.vSubsel(self.spikeCut.wfs(spikeAligned.alignIdx,:),...
                        nC, idxstart:idxstop);    

                    spikeAligned.tau = zeros(size(spikeAligned.unalignedwfs,1),1);
                    spikeAligned.maxIdx = self.P.spikeAlignment.cutLeft;

                    spikeAligned.restrictToIdx = spikeAligned.maxIdx-2:spikeAligned.maxIdx+5;
                    tic
                    if strcmp(self.P.spikeAlignment.method, 'onMax')
                        [spikeAligned.wfs spikeAligned.tau] = ...
                            mysort.wf.vAlignOnMax(spikeAligned.unalignedwfs, nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                    elseif strcmp(self.P.spikeAlignment.method, 'onUpsampledMean')
                        [spikeAligned.wfs spikeAligned.tau] = ...
                            mysort.wf.vAlignOnUpsampleMean(spikeAligned.unalignedwfs, nC,...
                            'maxIdx', self.P.spikeAlignment.maxIdx,...
                            'maxIter', self.P.spikeAlignment.maxIterations, ...
                            'initAlignment', self.P.spikeAlignment.initAlignment);   
                    elseif strcmp(self.P.spikeAlignment.method, 'onUpsampledMax')
                        [spikeAligned.tau spikeAligned.wfs] = ...
                            mysort.util.alignWaveformsUpsampleMax(spikeAligned.unalignedwfs, nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx, 'nIter', 2);
                    elseif strcmp(self.P.spikeAlignment.method, 'onUpsampledMin')
                        [spikeAligned.tau spikeAligned.wfs] = ...
                            mysort.wf.alignWaveformsUpsampleMax(-spikeAligned.unalignedwfs, nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx, 'nIter', 2); 
                        spikeAligned.wfs = -spikeAligned.unalignedwfs;
                    elseif strcmp(self.P.spikeAlignment.method, 'onMin')
                        [spikeAligned.wfs spikeAligned.tau] = ...
                            mysort.wf.vAlignOnMax(-spikeAligned.unalignedwfs, nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                        spikeAligned.wfs = -spikeAligned.wfs;
                    elseif strcmp(self.P.spikeAlignment.method, 'onAverageMax')
                        [spikeAligned.tau spikeAligned.wfs] = ...
                            mysort.wf.vAlignOnAverageMaxSample(spikeAligned.unalignedwfs, nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                    elseif strcmp(self.P.spikeAlignment.method, 'onAverageMin')
                        [spikeAligned.tau spikeAligned.wfs] = ...
                            mysort.wf.vAlignOnAverageMaxSample(-spikeAligned.unalignedwfs, nC,...
                            'maxIdx', spikeAligned.maxIdx, 'restrictToIdx', spikeAligned.restrictToIdx);
                        spikeAligned.wfs = -spikeAligned.wfs;
                    elseif strcmp(self.P.spikeAlignment.method, 'none')
                        spikeAligned.wfs = spikeAligned.unalignedwfs;
                    else
                        error(['unkown alignement method' self.P.spikeAlignment.method]);
                    end
                    toc
                    % cut away alignement artefacts at the ends
    %                 spikeAligned.unalignedwfs([1:max(1,max(spikeAligned.tau)) end+min(0,min(spikeAligned.tau))+1:end],:) = [];
                    save(self.files.spike_aligned_file, 'spikeAligned', '-v7.3');
                    disp('Done.')
                end
            end
            self.spikeAligned = spikeAligned; clear spikeAligned;        
        end

        %% --------------------------------------------------------------------
        function prewhitenSpikes(self)
            if exist(self.files.prewh_spike_file, 'file')
                disp('spikes already prewhitened...');
                spikePrewhitened = [];
                load(self.files.prewh_spike_file);
            else
                disp('prewhitening spikes...');
                % Cut the aligned waveforms to get rid of alignement artefacts
                % and reduce to the final waveform that will be used for
                % feature extraction
                idxstart = 1 + self.P.spikeAlignment.cutLeft - self.P.featureExtraction.cutLeft;
                idxstop  = idxstart + self.P.featureExtraction.Tf - 1;
                spikePrewhitened.wfs = mysort.wf.vSubsel(self.spikeAligned.wfs, nC, idxstart:idxstop);      
                % Build the noise covariance matrix and load it
                ccol_loaded = self.noise.Cestself.CCol;
                ccol_loaded(1:nC, 1:nC) = ccol_loaded(1:nC, 1:nC) + diag(diag(ccol_loaded(1:nC, 1:nC)));
                spikePrewhitened.ccol_loaded = ccol_loaded/2;
                spikePrewhitened.C = mysort.noise.ccol2Cte(spikePrewhitened.ccol_loaded, self.P.featureExtraction.Tf);
                spikePrewhitened.U = chol(spikePrewhitened.C);
                % Prewhiten
                spikePrewhitened.wfs = spikePrewhitened.wfs/spikePrewhitened.U; 
                save(self.files.prewh_spike_file, 'spikePrewhitened', '-v7.3');
                disp('Done.')    
            end
            self.spikePrewhitened = spikePrewhitened; clear spikePrewhitened;
        end

        %% --------------------------------------------------------------------
        function fetExtraction(self)
            if exist(self.files.fet_spike_file , 'file')
                disp('features already calculated...');
                spikeFeatures = [];
                load(self.files.fet_spike_file );
            else
                disp('Calculating features...');
                spikeFeatures.X = mysort.util.dimReductionPCA(self.spikePrewhitened.wfs,...
                    self.P.featureExtraction.nDims, [], 3*1000000);
                save(self.files.fet_spike_file, 'spikeFeatures', '-v7.3');
                disp('Done.')    
            end
            self.spikeFeatures = spikeFeatures; clear spikeFeatures;
        end

        %% --------------------------------------------------------------------
        function runMeanShift(self)
            if exist(self.files.meanshift_spike_file, 'file')
                disp('Already clustered...');
                clustering = [];
                load(self.files.meanshift_spike_file);
            else
                disp('Mean Shift Clustering...');
                if ~isempty(self.P.clustering.meanShiftBandWidth)
                    clustering.bandwidth = self.P.clustering.meanShiftBandWidth;
                else
                    clustering.bandwidth = sqrt(1.1*self.P.featureExtraction.nDims);
                end
                nSpFet = size(self.spikeFeatures.X,1);
                if self.P.clustering.maxSpikes < nSpFet
                    clusterIdx = randperm(nSpFet);
                    clustering.clusterIdx = sort(clusterIdx(1:self.P.clustering.maxSpikes));
                else
                    clustering.clusterIdx = 1:nSpFet;
                end
                X = self.spikeFeatures.X(clustering.clusterIdx,:)';
                tic
                [clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(X, clustering.bandwidth);
                toc  
                [clustering.ids clustering.clusterCenter] = MeanShiftClusterBundleResult(X', clustMembsCell, self.P.clustering.minSpikesPerCluster);
                clustering.classes = unique(clustering.ids);
                clusteredAlignedSpikes = self.spikeAligned.wfs(clustering.clusterIdx,:);
                clusteredCutSpikes = self.spikeCut.wfs(self.spikeAligned.alignIdx(clustering.clusterIdx),:);

                clustering.templatesCut = mysort.util.calculateClassMeans(clusteredCutSpikes, clustering.ids, 'usemedian');
                clustering.templatesAligned = mysort.util.calculateClassMeans(clusteredAlignedSpikes, clustering.ids, 'usemedian');
                save(self.files.meanshift_spike_file, 'clustering', '-v7.3');
            end
            self.clustering = clustering; clear clustering;
        end
        %% --------------------------------------------------------------------
        function runBOTMMatching(self)
            if exist(self.files.botm_matching_file, 'file')
                disp('botm matched with lda...');
                clusteringMatched = [];
                load(self.files.botm_matching_file);
            else
                disp('BOTM matching...');
                nSpCut = size(self.spikeCut.wfs,1);
                clusteringMatched.ids = zeros(nSpCut,1);
                clusteringMatched.template_ids = 0;
                clusteringMatched.spikeCutAligned = [];
                clusteringMatched.templates = [];
                clusteringMatched.maxTausPerSpike = zeros(nSpCut,1);
                clusteringMatched.ts = self.spikeDetectionMerged.ts(self.spikeCut.cutIdx);

                if size(self.clustering.templatesCut,1) < 2
                    disp('Cannot run BOTM Matching, not enough templates!');
                else
                    T = self.clustering.templatesAligned(2:end,:); % ignore noise template 
                    nT = size(T,1);
                    C = (self.noise.C_time_aligned +diag(diag(self.noise.C_time_aligned)))/2;

                    % Do matching only on temporal window defined by the aligned
                    % templates
                    idxstart = 1 + self.P.spikeCutting.cutLeft - self.P.spikeAlignment.cutLeft;
    %                 idxstop  = idxstart + self.P.spikeAlignment.Tf - 1;

                    [D maxTaus TupShiftDown FupShiftDown tauRange tauIdx subTauRange subTauIdx] = mysort.wf.vTemplateMatching(...
                        self.spikeCut.wfs, T, nC, idxstart, 'maxShift', 5, 'upsample', 5, 'noiseCovariance', C);
                    F = FupShiftDown(:,:,1);
                    EF = diag(T*F');                     % compute energies
                    Prior  = self.P.templateMatchingCut.prior;                
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

                    if 0 %~isempty(self.P.templateMatchingCut.residualPeakSmallerThanStdNoise)
                        % compute the "explained energy"
                        notNoiseIdx = find(notNoiseIdx);
                        nSmatched = length(notNoiseIdx);
                        E = zeros(size(self.spikeCut.wfs));
                        tTf = size(TupShiftDown,2);
                        for s = 1:8 %nSmatched
                            myId = clusteringMatched.ids(notNoiseIdx(s));
                            mySubTauIdx = clusteringMatched.maxSubTausIdxPerSpike(notNoiseIdx(s));
                            myTau = floor(clusteringMatched.maxTausPerSpike(notNoiseIdx(s)));
                            myStartIdx = idxstart-myTau;
                            E(notNoiseIdx(s),myStartIdx:myStartIdx+tTf-1) = squeeze(TupShiftDown(myId,:,mySubTauIdx));
                        end
                        R = self.spikeCut.wfs-E;
                        figure;
                        subplot(2,1,1);
                        plot(self.spikeCut.wfs(1:4,:)');
                        hold on
                        plot(E(1:4,:)');
                        subplot(2,1,2); plot(R(1:4,:)');
                    end

                    % shift the spikes to compute the templates
                    clusteringMatched.spikeCutAligned = mysort.wf.vShift(self.spikeCut.wfs, nC, -round(clusteringMatched.maxTausPerSpike), 1);
                    clusteringMatched.templates = mysort.util.calculateClassMeans(...
                        clusteringMatched.spikeCutAligned, clusteringMatched.ids, 'usemedian');  
                    clusteringMatched.templatesIDs = unique(clusteringMatched.ids);
                    clusteringMatched.ts = self.spikeDetectionMerged.ts(self.spikeCut.cutIdx)+clusteringMatched.maxTausPerSpike;
                    if 0
                        N = 1000;
                        figure; subplot(2,1,1); plot(self.spikeCut.wfs(1:N,:)');
                        subplot(2,1,2);
                        plot(clusteringMatched.spikeCutAligned(1:N,:)');
                    end
                end
                save(self.files.botm_matching_file, 'clusteringMatched', '-v7.3'); 
                disp('Done.')             
            end
            self.clusteringMatched = clusteringMatched;
        end
        %% --------------------------------------------------------------------
        function mergeClustersMS(self)
            if ~self.P.mergeTemplates.merge
                return
            end        
            if exist(self.files.merge_ms_clusters, 'file')
                disp('Already merged...');
                clusteringMerged = [];
                load(self.files.merge_ms_clusters);
            else   
                disp('Merging MS clusters...');

                % in self.clusteringMatched.ids the ids are going from 0 to nT with
                % 0 being not matched and noise, thus there is one "real" template less
                % than units (if any spike has an ID of 0 !!)
                realTemplateIdx = find(self.clusteringMatched.template_ids>0);
                nT = length(realTemplateIdx);
                if nT == 0
                    % we have only the noise cluster
                    D = [];
                    maxT = [];
                    groups = {};
                else
                    % resample only templates that do not have id==0
                    resampledTemplates = mysort.util.resampleTensor(mysort.wf.v2t(...
                        self.clusteringMatched.templates(realTemplateIdx,:), nC), 3, 1);
    %                 mysort.plot.waveforms(resampledTemplates, 'IDs', self.clusteringMatched.template_ids(realTemplateIdx), 'stacked', 0);

                    % the indices in the groups are indices into
                    % realTemplateIdx which is index into self.clusteringMatched.template_ids
                    [groups maxT D] = ana.mergeTemplates(resampledTemplates, self.noise.meanNoiseStd,...
                        'maxRelativeDistance', self.P.mergeTemplates.ifMaxRelDistSmallerPercent/100,...
                        'minCorrelation', self.P.mergeTemplates.atCorrelation);
                end
                clusteringMerged.D = D;
                clusteringMerged.maxT = maxT;
                clusteringMerged.groups = groups;            
                clusteringMerged.templates = zeros(length(groups), self.P.spikeCutting.Tf*nC);
                clusteringMerged.ids = zeros(length(self.clusteringMatched.ids),1);

                % we need to merge now the indices in ids accoring to the groups in
                % groups. But groups points into realTemplateIdx !!
                for i=1:length(groups)
                    group_ids = self.clusteringMatched.template_ids(realTemplateIdx(clusteringMerged.groups{i}));
                    idx = ismember(self.clusteringMatched.ids, group_ids);
                    clusteringMerged.ids(idx) = i;
                    clusteringMerged.templates(i,:) = median(self.clusteringMatched.spikeCutAligned(idx,:),1);
                end 
    %             mysort.plot.waveforms(self.spikeAligned.wfs, 'nC', nC, 'IDs', clusteringMerged.ids, 'stacked', 0)

                save(self.files.merge_ms_clusters, 'clusteringMerged', '-v7.3');
                fprintf('Templates after merging: %d\n', length(clusteringMerged.groups));
            end
            self.clusteringMerged = clusteringMerged;
        end    

        %% --------------------------------------------------------------------
        function runBOTM(self)
            if ~self.P.botm.run
                return
            end

            if exist(self.files.botm_file, 'file')
                disp('already sorted with botm...');
                botm = [];
                load(self.files.botm_file);
            else
                disp('BOTM sorting...');
                T = self.clusteringMerged.templates;
                Cest = self.noise.CestS;
                Cest.CCol(1:nC,:) = Cest.CCol(1:nC,:) + diag(diag(Cest.CCol(1:nC,:)));
                Cest.CCol = Cest.CCol/2;
                Cest = mysort.noise.Covest2(Cest);
                idxstart = 1 + self.P.spikeCutting.cutLeft - self.P.botm.cutLeft;
                idxstop  = idxstart + self.P.botm.Tf - 1;
                botm.templates = mysort.wf.vSubsel(T, nC, idxstart:idxstop);
                botmsorter = mysort.sorters.BOTM(Cest, self.P.botm.Tf, botm.templates,...
                            'upsample', 3, 'spikePrior', self.P.botm.prior,...
                            'max_num_of_channels_per_template', 30, 'adaptOnInit', 1);    

                botmsorter.DH = DS;    
                botm.gdfUnclean = botmsorter.sort(DS);
                % Remove artefacts from gdf
                botm.removeIdx = mysort.epoch.findPointsInEpochs(botm.gdfUnclean(:,2), self.artefactDetection.epochs)>0;               
                botm.gdf = botm.gdfUnclean(~botm.removeIdx,:);
                save(self.files.botm_file, 'botm');
                disp('Done.')    
            end
            self.botm = botm; 
        end
        %% --------------------------------------------------------------------
        function computeBOTMTemplates(self)
            if ~isempty(gdfbotm)
                if exist(self.files.botm_templates_file, 'file')
                    disp('loading botm templates...');
                    D = load(self.files.botm_templates_file);
                    templatesAfterBOTM = D.templatesAfterBOTM;
                    clear D
                else
                    disp('computing botm templates...');
                    uclasses = unique(gdfbotm(:,1));
                    templatesAfterBOTM = zeros(length(uclasses), self.P.spike_cut_Tf*nC);
                    for i=1:length(uclasses)
                        myts = gdfbotm(gdfbotm(:,1)==uclasses(i),2);
                        if size(myts,1) > 1000
                            myts = myts(1:1000);
                        end
                        spikes = self.DS.getWaveform(myts, self.P.spike_cut_cutLeft, self.P.spike_cut_Tf);
                        templatesAfterBOTM(i,:) = mean(spikes,1);
                    end
                    save(self.files.botm_templates_file, 'templatesAfterBOTM');
                    disp('Done.')    
                end            
            end
            self.botmDetails = botmDetails; clear botmDetails;
        end
        %% --------------------------------------------------------------------
        function residualArtefactDetection(self)
            if ~isempty(gdfbotm)
                if exist(self.files.resArtefact_file, 'file')
                    disp('residual artefacts already detected');
                    D = load(self.files.resArtefact_file);
                    gdfbotmcleaned = D.gdfbotmcleaned;
                    resArtefactEpochs = D.resArtefactEpochs; clear D;
                else
                    disp('Residual atrefact detection ...');
                    SpSoBOTM = mysort.spiketrain.SpikeSortingContainer('botm', gdfbotm, ...
                        'templateCutLeft', self.self.P.spike_cut_cutLeft, ...
                        'templateCutLength', self.self.P.spike_cut_Tf, ...
                        'templateWfs', mysort.wf.v2t(templatesAfterBOTM, nC),...
                        'wfDataSource', DS, 'nMaxSpikesForTemplateCalc', 1000);
                    self.DS.addSpikeSorting(SpSoBOTM); self.DS.setActiveSpikeSortingIdx('botm')            
                    self.DS.bReturnSortingResiduals = 1;
                    A = ana.moritzheimdahl.ResidualArtefactDetector(DS);
    %                 mysort.plot.SliderDataAxes({DS, A}, 'channelSpacers', [5000 0]);
                    resArtefactEpochs = A.getEpochs();
                    gdfbotmcleaned = gdfbotm;
                    removeIdx = mysort.epoch.findPointsInEpochs(gdfbotm(:,2), resArtefactEpochs)>0;
                    gdfbotmcleaned(removeIdx,:) = [];      
                    save(self.files.resArtefact_file, 'resArtefactEpochs', 'gdfbotmcleaned');
                    self.DS.bReturnSortingResiduals = 0;
                    disp('Done.')    
                end
            end
            self.botmResiduals = botmResiduals; clear botmResiduals;
        end
        %% --------------------------------------------------------------------
        function computeBOTMCleanedTemplates(self)
            if ~isempty(gdfbotmcleaned)
                if exist(self.files.botm_templates_cleaned_file, 'file')
                    disp('loading botm cleaned templates...');
                    D = load(self.files.botm_templates_cleaned_file);
                    templatesAfterBOTMcleaned = D.templatesAfterBOTMcleaned;
                    clear D
                else
                    disp('computing botm cleaned templates...');
                    uclasses = unique(gdfbotmcleaned(:,1));
                    templatesAfterBOTMcleaned = zeros(length(uclasses), self.P.spike_cut_Tf*nC);
                    for i=1:length(uclasses)
                        myts = gdfbotmcleaned(gdfbotmcleaned(:,1)==uclasses(i),2);
                        if size(myts,1) > 1000
                            myts = myts(1:1000);
                        end
                        spikes = self.DS.getWaveform(myts, self.P.spike_cut_cutLeft, self.P.spike_cut_Tf);
                        templatesAfterBOTMcleaned(i,:) = mean(spikes,1);
                    end
                    save(self.files.botm_templates_cleaned_file, 'templatesAfterBOTMcleaned');
                    disp('Done.')    
                end            
            end
            self.botmCleaned = botmCleaned; clear botmCleaned;                
        end
    end
end