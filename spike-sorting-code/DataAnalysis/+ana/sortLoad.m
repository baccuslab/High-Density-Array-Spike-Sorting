function [S P] = sortLoad(dpath, name, varargin)
%     P.dummy = 1;
%     P = mysort.util.parseInputs(P, varargin, 'error');
  
    % Define file names
    load([fullfile(dpath, name) '.P.mat'], 'S');
    P = S.P;
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
   
    disp('All loaded.!');
    

    %% --------------------------------------------------------------------
    function runArtefactDetection()
        artefactDetection.epochs = [];        
        if P.artefactDetection.use
            if exist(S.files.artefact_file, 'file')
                disp('artefacts already detected');
                load(S.files.artefact_file, 'artefactDetection');
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
        end
        S.spikeDetection = spikeDetection; clear spikeDetection;
    end

    %% --------------------------------------------------------------------
    function mergeDetectedSpikes()
        if exist(S.files.spike_det_merged_file, 'file')
            disp('spikes already merged');
            spikeDetectionMerged = [];
            load(S.files.spike_det_merged_file);
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
        end
        S.spikeCut = spikeCut; clear spikeCut;
    end

    %% --------------------------------------------------------------------
    function noiseEstimation()
        if exist(S.files.cov_file, 'file')
            disp('cov already estimated...');
            noise = [];
            load(S.files.cov_file);
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
        end
        S.spikePrewhitened = spikePrewhitened; clear spikePrewhitened;
    end

    %% --------------------------------------------------------------------
    function fetExtraction()
        if exist(S.files.fet_spike_file , 'file')
            disp('features already calculated...');
            spikeFeatures = [];
            load(S.files.fet_spike_file );
        end
        S.spikeFeatures = spikeFeatures; clear spikeFeatures;
    end
    
    %% --------------------------------------------------------------------
    function runMeanShift()
        if exist(S.files.meanshift_spike_file, 'file')
            disp('Already clustered...');
            clustering = [];
            load(S.files.meanshift_spike_file);
        end
        S.clustering = clustering; clear clustering;
    end
    %% --------------------------------------------------------------------
    function runBOTMMatching()
        if size(S.clustering.templates,1) < 2
            disp('Cannot run BOTM Matching, not enough templates!');
            return
        end
        if exist(S.files.botm_matching_file, 'file')
            disp('botm matched with lda...');
            clusteringMatched = [];
            load(S.files.botm_matching_file);
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
        end
        S.clusteringMerged = clusteringMerged;
    end    

    %% --------------------------------------------------------------------
    function runBOTM()
        S.botm = [];
        if ~P.botm.run
            return
        end    
        botm = [];
        if exist(S.files.botm_file, 'file')
            disp('already sorted with botm...');
            load(S.files.botm_file);
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
            end            
        end
        S.botmCleaned = botmCleaned; clear botmCleaned;                
    end
end