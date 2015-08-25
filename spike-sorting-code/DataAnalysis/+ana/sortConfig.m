function [rgcs gdf_merged T_merged localSorting localSortingID D] = sortConfig(dpath, hdmea, hdmea_session_idx, elGroupNumbers, elGroupIndices, runName, nCPUS, save2UMS, P)
% elgroups must contain the electrode numbers to be sorted together
if nargin < 7
    nCPUS = 2;
end
if nargin < 8
    save2UMS = 1;
end
if nargin < 9
    P = struct();
    P.spikeDetection.method = '+-';
    P.artefactDetection.use = 0;
    P.botm.run = 0;
    P.spikeCutting.maxSpikes = 100000000;
    
    P.spikeAlignment.initAlignment = '-';
    
    P.mergeTemplates.merge = 1;
    P.mergeTemplates.upsampleFactor = 3;
    P.mergeTemplates.atCorrelation = .92;
    P.mergeTemplates.ifMaxRelDistSmallerPercent = 30;
end
%% Run Local Sortings
tic;
if nCPUS > 1
    if matlabpool('size')==0
        disp('Opening Pool.')
        matlabpool(nCPUS);
    elseif matlabpool('size') ~= nCPUS
        disp('Re-Opening Pool.')
        matlabpool close;
        matlabpool(nCPUS);
    end
end

startOffsets = [];
if nCPUS > 1
    disp('Starting parfor')
    parfor i=1:length(elGroupIndices)
        elgroupNr = elGroupNumbers{i};
        elgroupIdx = elGroupIndices{i};
        fprintf('Starting %d\n', i);
        outpath = fullfile(dpath, sprintf('group%03d', i));
        if ~exist(outpath, 'file'); mkdir(outpath); end

        % make a local copy of the hdmea object, this is necessary for the
        % hd5 identifiers of another session parallelization since the hdf5 library cannot reuse the 
        hdmea_local = hdmea.copy();
        [S sortingStartOffsets] = ana.sortElectrodeGroup(hdmea_local, hdmea_session_idx, outpath, elgroupNr, elgroupIdx, runName, P);
        fprintf(['Done ' dpath ' group %d\n'], i);
        startOffsets(i) = sortingStartOffsets;
    end
elseif nCPUS == 1
    disp('Starting with normal for')
    for i=1:length(elGroupIndices)
        elgroupNr = elGroupNumbers{i};
        elgroupIdx = elGroupIndices{i};
        fprintf('Starting %d\n', i);
        outpath = fullfile(dpath, sprintf('group%03d', i));
        if ~exist(outpath, 'file'); mkdir(outpath); end

        % make a local copy of the hdmea object, this is necessary for the
        % hd5 identifiers of another session parallelization since the hdf5 library cannot reuse the 
        [S sortingStartOffsets] = ana.sortElectrodeGroup(hdmea, hdmea_session_idx, outpath, elgroupNr, elgroupIdx, runName, P);
        fprintf(['Done ' dpath ' group %d\n'], i);
        startOffsets(i) = sortingStartOffsets;
    end    
end
disp('+-+-+-+-+-+-+-+-+-+-+');
fprintf('Done with %s, run %s \n ', dpath, runName);
disp('+-+-+-+-+-+-+-+-+-+-+');
toc

%% Process all local sortings into a final sorting
[gdf_merged T_merged localSorting localSortingID] =...
    ana.processLocalSortings(dpath, runName, elGroupNumbers, elGroupIndices);
save(fullfile(dpath, [runName '_results.mat']), 'gdf_merged', 'T_merged', 'localSorting', 'localSortingID');


%% Build structure for export
if save2UMS
    disp('Getting Session Lengths...');
    sessionLengths =[];
    if isa(hdmea, 'mysort.ds.MultiSessionInterface');
        for i=1:length(hdmea_session_idx)
            sessionLengths(i) = size(hdmea.sessionList(hdmea_session_idx(i)),1);
        end    
        MEfull = hdmea.getAllSessionsMergedMultiElectrode();
    else
        sessionLengths = size(hdmea,1);
        MEfull = hdmea.MultiElectrode();
    end
    
    % mysort.plot.templates2D(T_merged, MEfull.electrodePositions, 3)    
    disp('Converting to UMS2000...');
    rgcs = ana.createUMS2000Structure(dpath, runName, hdmea,...
        hdmea_session_idx, sessionLengths, elGroupNumbers, MEfull, T_merged, ...
        localSorting, mod(localSortingID, 100));
    disp('Saving to UMS2000...');
    for i=1:length(rgcs)
        localSorting = rgcs{i};
        save(fullfile(dpath, [runName '_' sprintf('group%03d',i) '_Export4UMS2000.mat']), 'localSorting', 'startOffsets', '-v7.3')
    end
end
disp(['All done for this config. ' dpath]);

