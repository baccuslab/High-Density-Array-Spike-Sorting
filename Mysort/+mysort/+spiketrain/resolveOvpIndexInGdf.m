function [gdfresolved wasResolved]= resolveOvpIndexInGdf(gdf, ovpIdx)
    % gdf contains a gdf in which some units might represent overlapping
    % spikes of two other units. Those units are given in ovp index
    % ovpIdx is a matrix with 4 columns [gdfIDs unit1ID unit2ID shift]
    % the first column represents the unit id of the gds, i.e., the first
    % column of the gdf, which should be resolved into two units.
    % the second and third represent these two units.
    % the forth column is the shift of time of the second replacement unit.
    % Example:
    %   gdf = [1 100];
    %   ovpIdx = [1 2 3 -5];
    %   this means that spikes of unit 1 should be replaced by two spikes
    %   each, one of unit 2 with the same time point as the original spike
    %   of unit 1, and one spike of unit 3 with a shift of 5 samples
    %   leftwards
    %   gdfresolved = [3 95
    %                  2 100];
    if isempty(gdf) || isempty(ovpIdx)
        gdfresolved = [];
        wasResolved = [];
        return
    end
    
    gdfUnits = unique(gdf(:,1));
    newSpikes = [];
    gdfresolved = gdf;
    wasResolved = zeros(size(gdf,1),1);
    % check if any replacement unit should also be replaced. this cannot be
    assert(~any(ismember(ovpIdx(:,1), [ovpIdx(:,2); ovpIdx(:,3)])), 'A replacement unit should also be replaced, that cannot work!');
    for i = 1:size(ovpIdx(:,1))
        replacedUnit = ovpIdx(i,1); 
        newUnitWithoutShift = ovpIdx(i,2);
        newUnitWithShift    = ovpIdx(i,3);
        shift = ovpIdx(i,4);
        replacementSpikeIdx = find(gdf(:,1) == replacedUnit);
        nReplacedSpikes = length(replacementSpikeIdx);
        newSpikes = [newSpikes; [repmat(newUnitWithShift, nReplacedSpikes, 1) gdf(replacementSpikeIdx,2)-shift]];
        gdfresolved(replacementSpikeIdx,1) = newUnitWithoutShift;
        wasResolved(replacementSpikeIdx) = 1;
    end
    wasResolved = [wasResolved; ones(size(newSpikes,1),1)];
    [gdfresolved sortIdx] = sortrows([gdfresolved; newSpikes],2);
    wasResolved = wasResolved(sortIdx);