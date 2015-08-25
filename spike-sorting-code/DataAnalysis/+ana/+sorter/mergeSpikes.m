function allspikes = mergeSpikes(detUp, pksUp, detDown, pksDown, mergeSpikesMaxDist)
    detUp(cellfun(@isempty, pksUp)) = [];
    pksUp(cellfun(@isempty, pksUp)) = [];
    detDown(cellfun(@isempty, pksDown)) = [];
    pksDown(cellfun(@isempty, detUp)) = [];
    if isempty(detUp)
        allspikes = [cell2mat(detDown) cell2mat(pksDown)];        
    elseif isempty(detDown)
        allspikes = [cell2mat(detUp)   cell2mat(pksUp)];        
    else
        allspikes = [cell2mat(detUp)   cell2mat(pksUp);
                     cell2mat(detDown) cell2mat(pksDown)];
    end
    nSp = size(allspikes,1);
    
    if nSp == 0
        allspikes = [];
        return
    end
    
    % add the channel information to allspikes
    allspikes(nSp,3) = 0;
    nextIdx = 1;
    LL = max(length(detUp),...
             length(detDown));
    for i=1:LL
        nSpUp = 0;
        nSpDo = 0;
        if ~isempty(detUp)
            nSpUp = length(detUp{i});
        end
        if ~isempty(detDown)
            nSpDo = length(detDown{i});
        end
        allspikes(nextIdx:nextIdx+nSpUp+nSpDo-1,3) = i;
        nextIdx = nextIdx+nSpUp+nSpDo;
    end
    allspikes = sortrows(allspikes,1);
    fprintf('%d spikes found total before merging\n', nSp);

    allspikes  = mysort.spiketrain.mergeSingleElectrodeDetectedSpikes(allspikes, ...
        mergeSpikesMaxDist);
