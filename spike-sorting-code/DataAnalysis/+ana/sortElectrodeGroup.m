function [S sortingStartOffsets] = sortElectrodeGroup(hdmea, hdmea_sessionIdx, outpath, elGroupNumbers, elGroupIndices, runName, P)
    % this assumes that all sessions in sessionIdx share the same electrode
    % configuration !
    % todo, check for that!
    S = struct();
%     MAXDATALENGTH = 1000*1*60*20000;
%     OMMITBEGINNING = 0; %30*20000;
%     me = hdmea.getMergedMultiElectrode4Sessions(hdmea_sessionIdx);
%     elGroupIndices = me.getElIdx4ElNumber(elGroupNumbers);
    disp('Loading Data');
    nS = length(hdmea_sessionIdx);
    sessionLength = []; %hdmea.getSessionLength(hdmea_sessionIdx);
    sortingStartOffsets = zeros(1,nS);
%     if sum(sessionLength) > MAXDATALENGTH
%         % data is too long.
%         lps = floor(MAXDATALENGTH/nS);
%         X = [];
%         for i=1:nS
%             s1 = max(1,min(sessionLength(i)-lps, OMMITBEGINNING));
%             s2 = min(sessionLength(i), s1+lps);
%             X = [X; hdmea(s1:s2, elGroupIndices, hdmea_sessionIdx(i))];
%             sortingStartOffsets(i) = s1;
%         end
%     else
%         X = hdmea(:, elGroupIndices, hdmea_sessionIdx);
%     end
    if isa(hdmea, 'mysort.ds.MultiSessionInterface')
        DS = hdmea.getSessionObject(hdmea_sessionIdx);
    else
        DS = hdmea;
    end
    DS.restrictToChannels(elGroupIndices);
    if isa(DS, 'mysort.ds.PreProcessedDataSourceInterface')
        DS.bDisablePreprocessing = true;
    end
%     DS = mysort.ds.Matrix(X, hdmea.samplesPerSecond);

    %% Start Sorting
    S = ana.sort(DS, outpath, runName, P);
    
    %% Re-Estimate templates on all electrodes of that config after sorting
    % reset DS to be the DataSource with all channels
    DS.restrictToChannels([]);
    
    templateFile = fullfile(outpath, [runName '_templates.mat']);
%     if exist(templateFile, 'file')
%         delete(templateFile)
%     end
    if ~exist(templateFile, 'file')
        disp('Estimating Templates...');
        L = [];
        if isa(hdmea, 'mysort.ds.MultiSessionInterface');
            for i=1:length(hdmea_sessionIdx)
                L(i) = size(hdmea.sessionList(hdmea_sessionIdx(i)),1);
            end
        else
            L = size(hdmea,1);
        end
        if isempty(S.clusteringMerged.ids) || isempty(S.clusteringMatched.ts)
            gdf = [];
        else
            gdf = [S.clusteringMerged.ids S.clusteringMatched.ts(:)];
        end
    %     gdf = mysort.spiketrain.gdf2multiSessionGdf(S.botm.gdf(S.botm.gdf(:,1)>0,:), [0 cumsum(L)], hdmea_sessionIdx);
    %     bwfmngr = mysort.wf.MultiSessionBufferedWfManager(hdmea, gdf(:,2), gdf(:,1), gdf(:,3), 10, 55, 1:size(gdf,1), hdmea.getAllSessionsMergedMultiElectrode);
    %     mysort.mea.Template
        sortingContainerFile = fullfile(outpath, [runName '_SortingContainer.h5']);
        if exist(sortingContainerFile, 'file')
            delete(sortingContainerFile);
        end
        if isempty(gdf)
            gdfList = cell(1,length(hdmea_sessionIdx));
        else
            gdfList = mysort.spiketrain.splitGdf(gdf(gdf(:,1)>0,:), [0 cumsum(L)], hdmea_sessionIdx);
            for i=1:length(gdfList)
                if ~isempty(gdfList{i})
                    gdfList{i}(:,2) = gdfList{i}(:,2)+sortingStartOffsets(i);
                end
            end
        end
        if isa(hdmea, 'mysort.ds.MultiSessionInterface');
            MSSpContainer = mysort.mea.PersistentMultiSessionSpikeSortingContainer(sortingContainerFile, 'mssort', gdfList);
            MSSpContainer.computeTemplates(hdmea);
            MSSpContainer.save();
            ME = hdmea.getAllSessionsMergedMultiElectrode();
            nSourceSpikesPerTemplateAndChannel = MSSpContainer.TemplateManager.getNSourceSpikes4MultiElectrode(ME);
            [wfs cutleft] = MSSpContainer.TemplateManager.getWaveforms4MultiElectrode(ME);

        else
            MSSpContainer = mysort.spiketrain.PersistentSpikeSortingContainer(sortingContainerFile, 'bla', 'mssort', gdfList{1}, 'wfDataSource', hdmea);
            MSSpContainer.save2File();
            MSSpContainer.computeTemplates(hdmea);
            ME = hdmea.MultiElectrode();
            nSourceSpikesPerTemplateAndChannel = [];
            wfs = MSSpContainer.templateWfs;
            cutleft = MSSpContainer.templateCutLeft;
        end
        
        
        MES = ME.toStruct();
        save(templateFile, 'wfs', 'cutleft', 'MES', 'elGroupIndices', 'nSourceSpikesPerTemplateAndChannel'); 
    else
        disp('Templates already estimated');
    end
    
%     hdmea.addSpikeSorting(MSSpContainer);
    
%     mysort.plot.SliderDataAxes(hdmea, 'channelSpacers', 100);
%     gdf51 = S.botm.gdf(ismember(gdf(:,1), [1 5]),:);
%     mysort.plot.isiViolations([], [ones(size(gdf51,1),1) gdf51(:,2)], 20)
%     mysort.plot.xcorr(S.botm.gdf(ismember(S.botm.gdf(:,1), [1 5]),:), 'srate', 20000, 'binSize', 1);
