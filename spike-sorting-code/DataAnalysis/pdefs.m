function def = pdefs()
    bIsOffline = true; %false; %true;
%     this_path = pwd;
%     offlineFile = fullfile(this_path, 'offline_status_indicator.txt');
%     if exist(offlineFile, 'file')
%         bIsOffline = true;
%     end
    def.pdefsPath = pwd;
%     callerName = mysort.util.getCallerName();
    def.callerPath = [];
    if ~isempty(strfind(computer, 'WIN'))
        def.sortingOutPath = 'C:\localData\sortingOut\';
        def.networkTempShare = 'S:\group\hierlemann\Temp\FelixFranke\';
        def.localData = 'C:\localData\';
        def.localDataOnNetwork = fullfile(def.networkTempShare, 'LocalData');
        def.serverData = 'S:\group\hierlemann\recordings\collaborations\';
        def.serverDataRoska = 'S:\group\hierlemann\recordings\HiDens\Roska';  
        def.svnDataRoska = '';
        def.ravaDropBox = 'C:\Users\frankef\Dropbox\Paper NoiseCorrelations\';
        def.ravaSimulations = fullfile(def.localDataOnNetwork, 'RavaSimulations');
        def.serverSortingOut = 'S:\group\hierlemann\recordings\SpikeSortingOut\';
        def.micheleDSPaperDataPath = 'S:\group\hierlemann\recordings\HiDens\Roska\Fiscella_2014\backup_results';
        if bIsOffline
            def.ravaSimulations = 'C:\LocalData\RavaSimulations\';
            def.localDataOnNetwork = def.localData;
            def.micheleDSPaperDataPath = 'C:\LocalData\Michele\backup_results';
        end
        
    elseif ~isempty(strfind(computer, 'MACI64'))
        %def.sortingOutPath = '/Volumes/hierlemann/recordings/HiDens/SpikeSorting/';
        def.networkTempShare = '/Volumes/hierlemann/Temp/FelixFranke/';
        
        def.localData = '/Users/rolandd/tmp/';
        def.serverData = '/Volumes/hierlemann/recordings/collaborations/';        
        def.localDataOnNetwork = fullfile(def.networkTempShare, 'LocalData');
        %def.hamsterDS = fullfile('/Volumes', 'hierlemann', 'recordings', 'Mea1k', 'shared', '140821');
        
        def.Mea1kShared = fullfile('/Volumes', 'hierlemann', 'recordings', 'Mea1k', 'shared');
        def.mea1kRoland = fullfile('/Volumes', 'hierlemann', 'recordings', 'Mea1k', 'rolandd');
        
        def.recordingsRoland =  fullfile('/Volumes', 'rolandd$', 'hima-storage', 'recordings', 'Mea1k', 'rolandd');  %fullfile('/Volumes', 'hierlemann', 'recordings', 'Mea1k', 'rolandd'); % to be changed
        def.intermediateRoland =  fullfile('/Volumes', 'rolandd$', 'hima-storage', 'intermediate_data', 'Mea1k', 'rolandd'); % fullfile('/Volumes', 'hierlemann', 'intermediate_data', 'Mea1k', 'rolandd'); % to be changed
        def.analysedRoland = fullfile('/Volumes', 'hierlemann', 'AnalyzedData', 'Mea1k', 'rolandd');
        
        %def.serverDataRoska = '/links/groups/hima/recordings/HiDens/Roska/';
        %def.svnDataRoska = '/home/frankef/bel.svn/hima_internal/cmosmea_recordings/trunk/Roska';        
        %def.ravaDropBox = fullfile(def.localData, 'RavaDropboxLinux');
        %def.ravaSimulations = fullfile(def.localData, 'RavaSimulations');
        %def.serverSortingOut = '/net/bs-filesvr01/export/group/hierlemann/recordings/SpikeSortingOut/';
        %def.micheleDSPaperDataPath = '/links/groups/hima/recordings/HiDens/Roska/Fiscella_2014/backup_results';  
    else
        def.sortingOutPath = '/links/groups/hima/recordings/HiDens/SpikeSorting/';
        def.networkTempShare = '/links/groups/hierlemann/Temp/FelixFranke/';
        def.localData = '/net/bs-filesvr01/export/group/hierlemann/Temp/FelixFranke/LocalData/';
        def.localDataOnNetwork = fullfile(def.networkTempShare, 'LocalData');
        def.serverData = '/links/groups/hima/recordings/collaborations/';        
        def.serverDataRoska = '/links/groups/hima/recordings/HiDens/Roska/';
        def.mea1kData = '/links/groups/hima/recordings/Mea1k/';
        def.mea1kIntermediate = '/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/';
        
%        def.mea1kRoland = '/links/groups/hima/recordings/Mea1k/rolandd/';
        def.recordingsRoland = '/net/bs-filesvr02/export/group/hierlemann/recordings/Mea1k/rolandd/';%'/links/groups/hima/recordings/Mea1k/rolandd/'; % to be changed
        def.intermediateRoland = '/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/rolandd/'; %'/links/groups/hima/intermediate_data/Mea1k/rolandd/'; % to be changed
        def.analysedRoland = '/links/groups/hima/AnalyzedData/Mea1k/rolandd';
        
        def.svnDataRoska = '/home/frankef/bel.svn/hima_internal/cmosmea_recordings/trunk/Roska';        
        def.ravaDropBox = fullfile(def.localData, 'RavaDropboxLinux');
        def.ravaSimulations = fullfile(def.localData, 'RavaSimulations');
        def.serverSortingOut = '/net/bs-filesvr01/export/group/hierlemann/recordings/SpikeSortingOut/';
        def.micheleDSPaperDataPath = '/links/groups/hima/recordings/HiDens/Roska/Fiscella_2014/backup_results';
        %def.hamsterDS = fullfile('/home', 'rolandd', 'hima', 'recordings', 'Mea1k', 'shared', '140821');
        def.Mea1kShared = fullfile('/home', 'rolandd', 'hima', 'recordings', 'Mea1k', 'shared');
        
    end