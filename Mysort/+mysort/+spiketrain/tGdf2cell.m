
function [spike_trains tIDs uIDs] = tGdf2cell(tGDF)
    spike_trains = {};
    if isempty(tGDF)
        return
    end
    tIDs = unique(tGDF(:,1));
    uIDs  = unique(tGDF(:,2));
    nT = length(tIDs);
    nU = length(uIDs);
    spike_trains = cell(nT, nU);
    for trial = 1:nT
        for unit = 1:nU
            spike_trains{trial, unit} = tGDF(tGDF(:,1)==tIDs(trial) & tGDF(:,2)==uIDs(unit),3);
        end
    end