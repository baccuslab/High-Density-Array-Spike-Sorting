function eT = alignGetErrorTypes(R, sgdf, targetUnitID, unitIDs)
    % This function computes the error types from the alignment result
    % structure R for one unit "targetUnitID". The other neurons are all
    % treated as TP, except they have clustering errors. This is good if 
    % only one neurons ground truth is known and its errors should be 
    % plotted. sgdf is the sorting result gdf
    %
    % CAREFUL, the first ID in R.spikeLabel2 is considered to be 0, not 1!
    % if that is not the case use unitIDs to specify the unitsIDs for the
    % corresponding cells in R.spikeLabel2.
    if nargin < 4 || isempty(unitIDs)
        unitIDs = 0:length(R.spikeLabel2)-1;
    end
    def = mysort.util.defs();
    nSP = sum(R.nSP2);
    eT = def.tp*ones(nSP,1);
    for uidx = 1:length(unitIDs)
        i = unitIDs(uidx);
        eT(sgdf(:,1)==i) = R.spikeLabel2{uidx};
        if i ~= targetUnitID
            % This is not the target unit. Set all FPs to TP
            idx = find(sgdf(:,1)==i);
            idx = idx(R.spikeLabel2{uidx}==def.fp);
            eT(idx) = def.tp;
        end
    end