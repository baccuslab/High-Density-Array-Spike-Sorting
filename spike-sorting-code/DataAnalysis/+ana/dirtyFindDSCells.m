% function dirtyFindDSCells(ntk_filelist, outputPath)
    
    cd('/home/frankef/bel.svn/hima_internal/cmosmea_recordings/trunk/Roska/19Apr2012_DSGCs_rabbit/Matlab');
    outputPath = '/home/frankef/bel.svn/hima_internal/cmosmea_recordings/trunk/Roska/19Apr2012_DSGCs_rabbit/Matlab/DirtySortings';
    confs = load('configurations.mat');
    c = confs.configs(2);
    ntk_filelist = confs.flist(c.flistidx);
    system(['ls -lah ../proc/' ntk_filelist{1}])
    
    % This function does a "supposedly" quick spike sorting of the ntk
    % files given in ntk_filelist and uses hardcoded stimulus information
    % to estimate the direction selectivity index of the neurons found. 
    % Results are stored in outputPath

    cputimes = []; ti=1; tstart = tic;
    
    fprintf('Init...\n')
    sps = 20000; % samples per second
    hpf = 300;
    lpf = 7000;
    forder = 4;
    name = 'dirtySort';
    runName = 'r2';
    
    cutleft = 10;
    Tf = 45;
    
    nFiles = length(ntk_filelist);
    bufferLengthPerFileInMin = 1;
    bufferLengthPerFileInSamples = bufferLengthPerFileInMin*60*sps;
    % The individual files should not contain much more samples then estimated
    % in this way!
    
    % get number of channels
    nSamplesPerFiles = zeros(1, nFiles);   
    
    % Load the whole data of all ntk files
    siz = bufferLengthPerFileInSamples*2;
    fprintf('Loading raw data...\n')
    for i=1:2 %nFiles
        ntk=initialize_ntkstruct(ntk_filelist{i}, 'nofilters');
        [high_density_data ntk] = ntk_load(ntk, siz, 'images_v1');
        nC = size(high_density_data.sig, 2);
        if i==1
            % init data matrix
            X = zeros(nFiles*bufferLengthPerFileInSamples, nC);
            FR = zeros(nFiles*bufferLengthPerFileInSamples, 1);
        end
        
        nSamplesPerFiles(i) = size(high_density_data.sig,1);
        assert(nSamplesPerFiles(i) < siz, 'The file was longer than allowed! Only use this script for short files that contain the DS cell finding stimuli!');
        X (sum(nSamplesPerFiles(1:i-1))+1:sum(nSamplesPerFiles(1:i)),:) = high_density_data.sig;
        FR(sum(nSamplesPerFiles(1:i-1))+1:sum(nSamplesPerFiles(1:i)),1) = high_density_data.images.frameno;
    end
    X = X(1:sum(nSamplesPerFiles),:);
    FR = FR(1:sum(nSamplesPerFiles),:);
    stim_epochs = ana.michele.epochsFromFrameNo(FR);
    
%     figure;
%     plot(FR);
%     hold on
%     plot(stim_epochs(:,1), FR(stim_epochs(:,1)), 'dg');
%     plot(stim_epochs(:,2), FR(stim_epochs(:,2)), 'dr');
%     
    clear FR
    % filter the data and build DataSourceObject
    cputimes(ti) = toc(tstart); ti=ti+1; tstart = tic;
    fprintf('Prefiltering data...\n')
    hd = mysort.mea.filter_design(hpf, lpf, sps, forder);
    X = filtfilthd(hd, X);
    electrodePositions = [high_density_data.x' high_density_data.y']; 
    electrodeNumbers = high_density_data.channel_nr';
    
    % Make groupings of electrodes to be sorted independently
    [groupsidx nGroupsPerElectrode] = mysort.mea.constructLocalElectrodeGroups(electrodePositions(:,1), electrodePositions(:,2));
    % replace electrode indices with electrode numbers
    groups = {};
    for ii=1:length(groupsidx)
        groups{ii} = electrodeNumbers(groupsidx{ii});
    end    
    
    % sort the groups
    cputimes(ti) = toc(tstart); ti=ti+1; tstart = tic;
    P = struct();
    P.spikeDetection.method = '-';
    P.spikeDetection.thr = 4.5;
    P.artefactDetection.use = 0;
    P.botm.run = 0;
    P.spikeCutting.maxSpikes = 100000; % This is crucial for computation time
    
    P.noiseEstimation.minLength = 200000; % This also effects computation time
    P.noiseEstimation.minDistFromSpikes = 80;
    
    P.spikeAlignment.initAlignment = '-';   
    P.spikeAlignment.maxSpikes = 30000;
    P.clustering.maxSpikes = 30000;
    
    P.mergeTemplates.merge = 1;
    P.mergeTemplates.upsampleFactor = 3;
    P.mergeTemplates.atCorrelation = .92;
    P.mergeTemplates.ifMaxRelDistSmallerPercent = 30;
    % Prepare data
    XX = {};
    for ii=1:length(groupsidx)
        XX{ii} = X(:, groupsidx{ii});
    end
    DSFull = mysort.ds.Matrix(X, sps, name, electrodePositions, electrodeNumbers);    
    cputimes(ti) = toc(tstart); ti=ti+1; tstart = tic;
    parfor ii=1:length(groupsidx)
        elIdx = groupsidx{ii};
        DS = mysort.ds.Matrix(XX{ii}, sps, name, electrodePositions(elIdx,:), electrodeNumbers(elIdx));        
        [S P_] = ana.sort(DS, fullfile(outputPath, ['group' sprintf('%03d', ii)]), runName, P);
    end
    clear XX DS S;
    
    % Re-Estimate templates on all electrodes of that config after sorting
    cputimes(ti) = toc(tstart); ti=ti+1; tstart = tic;
    disp('Estimating Templates...');
    MES = DSFull.MultiElectrode.toStruct();
    for ii=1:length(groupsidx)
        matchingFile  = fullfile(fullfile(outputPath, ['group' sprintf('%03d', ii)]), [runName '.100botm_matching.mat']);
        sortingFile  = fullfile(fullfile(outputPath, ['group' sprintf('%03d', ii)]), [runName '.110clusters_meanshift_merged.mat']);
        templateFile = fullfile(fullfile(outputPath, ['group' sprintf('%03d', ii)]), [runName '_templates.mat']);
        S = load(sortingFile);
        M = load(matchingFile);
        units = unique(S.clusteringMerged.ids);
        wfs = zeros(Tf, size(DSFull,2), length(units));
        nSourceSpikesPerTemplateAndChannel = [];
        for uidx = 1:length(units)
            fprintf('.');
            tsIdx = find(S.clusteringMerged.ids == units(uidx));
            tsIdx = tsIdx(1:min(100, end));
            wfs_ = DSFull.getWaveform(round(M.clusteringMatched.ts(tsIdx)), cutleft, Tf);
            wfs(:,:,uidx) = mysort.wf.v2t(median(wfs_,1), size(DSFull,2));
            nSourceSpikesPerTemplateAndChannel(uidx, 1:size(DSFull,2)) = size(wfs_,1);
        end
        elGroupIndices = groupsidx{ii};
        save(templateFile, 'wfs', 'cutleft', 'MES', 'elGroupIndices', 'nSourceSpikesPerTemplateAndChannel'); 
        fprintf('.\n');
    end
    
    
%% Process all local sortings into a final sorting
cputimes(ti) = toc(tstart); ti=ti+1; tstart = tic;
disp('Postprocessing...');
[gdf_merged T_merged localSorting localSortingID] =...
    ana.processLocalSortings(outputPath, runName, groups, groupsidx);
save(fullfile(outputPath, [runName '_results.mat']), 'gdf_merged', 'T_merged', 'localSorting', 'localSortingID');

cputimes(ti) = toc(tstart); 
fprintf('Loading: %.3f\nFiltering: %.3f\nDataInit: %.3f\nSorting: %.3f\nEst Templates: %.3f\nPostproc: %.3f\n', cputimes);

%% Estimate direction selectivity
neurons = unique(gdf_merged(:,1));
nU = length(neurons);
nE = size(stim_epochs,1);
countsraw = zeros(nU, nE);
rates = counts;
for i=1:nU
    u = neurons(i);
    tsu = gdf_merged(gdf_merged(:,1)==u, 2);
    for e=1:nE
        countsraw(i,e) = sum(tsu>=stim_epochs(e,1) & tsu<=stim_epochs(e,2));
        rates(i,e) = countsraw(i,e)*sps/diff(stim_epochs(e,:));
    end
end


%% Sort counts 
[sortedDirs sortDirections] = sort(c.DIRECTIONS(1,:));
angvecs1 = [cos(deg2rad(sortedDirs))' sin(deg2rad(sortedDirs))'];
counts = countsraw(:,[sortDirections 8+sortDirections]);
rates  = rates(:,[sortDirections 8+sortDirections]);
normCounts = counts./repmat(sum(counts,2),1,nE);

%%
figure
axes
hold on

for i = 1:nU;
    u = neurons(i);
    tsu = gdf_merged(gdf_merged(:,1)==u, 2);
    plot(tsu, i*ones(1,length(tsu)), '.');
end
plot(stim_epochs(:,1), zeros(1,nE), 'dg');
plot(stim_epochs(:,2), zeros(1,nE), 'dr');

%%
dsidx = zeros(nU,1);
angles = dsidx;
vecsum = dsidx;
vectors1 = zeros(nU,2);
for i=1:nU
    % 1 2 3 4 5 6 7 8
    % |       |
    %   |       |
    %     |       |
    %       |       |
    cu = counts(i,1:8)+counts(i, 9:16);
    [maxc prefidx] = max(cu);
    if prefidx >= 5
        nullidx = prefidx-4;
    else
        nullidx = prefidx+4;
    end
    dsidx(i) = ( cu(prefidx)-cu(nullidx) )/( cu(prefidx)+cu(nullidx) );
    
    cun = cu/sum(cu);
    vectors1(i,:) = cun*angvecs1;
    angles(i) = atan2(vectors1(i,1), vectors1(i,2));
    vecsum(i) = norm(vectors1(i,1));
end

figure
bar([dsidx vecsum])
dsCells = find(dsidx>.3 & vecsum >.2);

%%
figure
axes
hold on
for i=1:length(dsCells)
    plot([0 vectors1(dsCells(i),1)], [0 vectors1(dsCells(i),2)],'k-');
end
    
%%
figure
axes; hold on
for i = 1:length(dsCells)
    u = neurons(dsCells(i));
    tsu = gdf_merged(gdf_merged(:,1)==u, 2);
    plot(tsu, i*ones(1,length(tsu)), '.');
end
plot(stim_epochs(:,1), zeros(1,nE), 'dg');
plot(stim_epochs(:,2), zeros(1,nE), 'dr');


%%
figure; 
subplot(2,1,1)
plot((rates(dsCells,1:8)+rates(dsCells, 9:16))'/2)
subplot(2,1,2)
plot((normCounts(dsCells,1:8)+normCounts(dsCells, 9:16))');


