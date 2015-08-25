
function [spikeTrains IDs] = checkSpikeTrainsOrGDF(st_gdf, P)
    IDs = [];
    spikeTrains = st_gdf;
    
    if isnumeric(st_gdf)
        if isempty(st_gdf)
            warning('no spike trains supplied. not plotting anything');
            return
        end
        [spikeTrains IDs] = mysort.spiketrain.fromGdf(st_gdf);
    end
    
    if isempty(P.IDs) && isempty(IDs)
        IDs = 1:size(spikeTrains,2);
    elseif ~isempty(P.IDs)
        IDs = P.IDs;
    end    