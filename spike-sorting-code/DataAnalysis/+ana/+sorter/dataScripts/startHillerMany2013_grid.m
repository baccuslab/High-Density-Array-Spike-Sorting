% matlabpool(12);
% clear all


% DONE:
% expNames = {'configs_2Feb2015_E_mea1k'};    % < 2 >
% expNames = {'configs_5Feb2015_A_mea1k_hamster_left_eye_left_periphery'
%             'configs_5Feb2015_B_mea1k_hamster_left_eye_bottom_periphery'
%             'configs_5Feb2015_C_mea1k_hamster_right_eye_left_periphery'
%             'configs_5Feb2015_D_mea1k_hamster_right_eye_bottom_periphery'}; % < 2 >
% expNames = {'configs_2Feb2015_C_mea1k'
%             'configs_2Feb2015_D_mea1k'
% expNames = {'configs_2Feb2015_A_mea1k'
%             'configs_2Feb2015_B_mea1k'
%             'configs_3Feb2015_A_mea1k'
%             'configs_3Feb2015_B_mea1k'
%             'configs_3Feb2015_C_mea1k'};    % < Postproc 2 >
% expNames = {'configs_29Jan2015_A_mea1k'
%             'configs_29Jan2015_B_mea1k'
%             'configs_30Jan2015_A_mea1k'
%             'configs_30Jan2015_B_mea1k'};     % < Postproc 4 > 
%             'configs_18Feb2014_mea1k'
%             'configs_4Feb2015_B_mea1k'
%             'configs_4Feb2015_C_mea1k'
%             'configs_4Feb2015_D_mea1k'
%             'configs_4Feb2015_E_mea1k'
%             'configs_4Feb2015_F_mea1k'
%             'configs_5Feb2015_A_mea1k'
%             'configs_5Feb2015_B_mea1k'            
% };    % < 0 >
% expNames = {'configs_11Feb2015_A_mea1k'
%             'configs_11Feb2015_B_mea1k'
%             'configs_11Feb2015_C_mea1k'}; 

% expNames = {'configs_rabbit_highdens'
%             'configs_rabbit_sparse'}; 
% expNames = {'configs_13Mar2015_A_mea1k'};  
% expNames = {'configs_17Mar2015_A_mea1k'};  
% expNames = {'configs_19Mar2015_DZ_mea1k'}; 
% 'configs_19Mar2015_A_mea1k' 
% 'configs_19Mar2015_B_mea1k'
% expNames = {'configs_19Mar2015_C_mea1k'}; % < 0 >
% expNames = {'configs_rabbit_highdens'};
% expNames = {'configs_rabbit_sparse'};
% expNames = {'configs_15Dec2014_B_bug_search'}; % < 2 >
% expNames = {'configs_bug_search_5Feb2015_B_mea1k_hamster_left_eye_bottom_periphery'};
% expNames = {'configs_17Mar2015_A_mea1k_rabbitV2'}; % <3>
% expNames = {'configs_19Mar2015_B_mea1k_rabbitV2'}; % <2>
% expNames = {'configs_19Mar2015_A_mea1k_rabbitV2'}; % <4>
% expNames = {'configs_16Apr2015_primate_mea1k'}; % <4>
% expNames = {'configs_Auto_vs_Manual_sorting'}; % <4>
% expNames = {'configs_19Mar2015_A_mea1k_rabbitV2'};
% expNames = {'configs_23Apr2015_primate_mea1k'};
% expNames = {'JAN_131030_mostActiveEls'}; h5Path = '/Sessions/000/sig';
% expNames = {'configs_25Jun2015_hamster_felix_A'
%             'configs_25Jun2015_hamster_felix_B'};
% expNames = {'configs_24Jun2015_hamster_felix'};
% expNames = {'configs_30Jun2015_daniel_133'
%             'configs_30Jun2015_daniel_134'};
% expNames = {'configs_7Jul2015_CNO_michele'};
% expNames = {'configs_8Jul2015_PSEM_left_michele'
%             'configs_8Jul2015_PSEM_right_michele'};
% expNames = {'configs_10Jul2015_CNO_michele'};
%             'configs_23Jul2015_daniel_141'
%             'configs_23Jul2015_daniel_142'
%             'configs_23Jul2015_daniel_143'
%             'configs_23Jul2015_daniel_146'
%             'configs_21Jul2015_hamster_felix_periphery'

% IN PROCESS:

expNames = {'configs_22Jul2015_hamster_felix_periphery'
            'configs_21Jul2015_hamster_felix_center'
            'configs_22Jul2015_hamster_felix_center'};
h5Path = [];

% QUEUE:
% expNames = {''};   


%%
pd = pdefs();

CONFIG_FILE_LOCATION = fullfile(pd.networkTempShare, 'Hillier_2013');
INTERMEDIATE_DATA_PATH = fullfile(pd.mea1kIntermediate, 'frankef');
submitForReal = true;
% submitForReal = false
runLocally = false;
bCorrectForMissingFrames = false;

%%
gridSpikeSorting();   


%% EXPORT FOR MICHELE TO SINGLE FOLDER

for eidx = 1:length(expNames)
    expName = expNames{eidx};
    pathdefs = PATHDEFS(eidx);
    if ~exist(pathdefs.exportResultPath, 'file')
        mkdir(pathdefs.exportResultPath)
    end
    MGDFS = {};
    R = {};
    DIFFFRAMES = {};
    for i=1:nConfigs
        sourceFile = fullfile(pathdefs.sortOutPaths{i}, [pathdefs.jobName '_results.mat'])
        assert(exist(sourceFile, 'file')>0, 'Source File not found!');
        D = load(sourceFile);

        % get individual files' length in samples
        FRAMES = {};
        
        if bIsMea1Krecording && bCorrectForMissingFrames
            for k=1:length(pathdefs.configNTKLists{i})
                F = mea1k.file(pathdefs.configNTKLists{i}{k});
                FRAMES{k} = F.getFrameNo();
                df = diff(FRAMES{k})-1;
                idx = find(df);
                DIFFFRAMES{k,i,1} = idx;
                DIFFFRAMES{k,i,2} = df(idx);
            end
        end
%          disp(['Align the sorting output with the original frame numbers. ' num2str(i) ' out of ' num2str(length(rawFiles))] )
    
        
        compoundMea = mysort.mea.compoundMea(pathdefs.h5FileLists{i}, 'useFilter', 0, 'name', 'PREFILT');
        L = compoundMea.X.getAllSessionsLength();

        % break gdf apart
        start_end_times = [0 cumsum(L)];
        assert(max(D.gdf_merged(:,2)) <= start_end_times(end), 'Spike Times out of Range !!');
        mgdf = mysort.spiketrain.gdf2multiSessionGdf(D.gdf_merged, start_end_times);
        MGDFS{i} = mgdf;
        for k=1:length(L)
            gdf = mgdf(mgdf(:,3)==k,[1:2]);
            if ~isempty(FRAMES)
                gdf(:,2) = FRAMES{k}(gdf(:,2));
            end
            R{k, i} = gdf;
        end
    end
    save(fullfile(pathdefs.inputPath, [expName '_resultsForMichele.mat']), 'MGDFS', 'pathdefs', 'R', 'bCorrectForMissingFrames', 'DIFFFRAMES');
end


return
% 
% if bIsMea1Krecording && bCorrectForMissingFrames
%     %% summary of missing frame numbers
%     N = cellfun(@(x) length(x), DIFFFRAMES(:,:,1));
%     fprintf('Number of blocks with missing frames: %d\n', sum(N));
%     NN = sum(cellfun(@(x) sum(x), DIFFFRAMES(:,:,2)));
%     fprintf('Number of missing frames total: %d\n', NN);
% else
%     disp('Frame number correction was not run!');
% end
