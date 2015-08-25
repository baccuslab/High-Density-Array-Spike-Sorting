function [xc, edges, binnedST, lags, xc_idx] = xcorrBinned(spikeTrains_or_gdf, binWidth, varargin)
    P.maxLag = 100;
    P.overwriteNumberOfUnitAssertion = false;
    P = mysort.util.parseInputs(P, varargin, 'error');
    
    if ~iscell(spikeTrains_or_gdf)
        st = mysort.spiketrain.fromGdf(spikeTrains_or_gdf);
        gdf = spikeTrains_or_gdf;
    else
        st = spikeTrains_or_gdf;
        gdf = mysort.spiketrain.toGdf(spikeTrains_or_gdf);
    end
    
    nUnits = length(st);
    if ~P.overwriteNumberOfUnitAssertion
        assert(nUnits < 200, 'Processing more than 200 spike trains is not recommended');
    end
    
    minTs = min(gdf(:,2));
    maxTs = max(gdf(:,2));
    edges = minTs:binWidth:(maxTs+binWidth) - binWidth/2;
    
    binnedST = mysort.spiketrain.toBins(st, edges);

    if isempty(P.maxLag)
        [xc, lags] = xcorr(binnedST', 'coeff');
    else
        [xc, lags] = xcorr(binnedST', P.maxLag, 'coeff');
    end
    
    xc_idx = zeros(2, size(xc,2));
    count = 0;
    for i=1:nUnits
        for j=1:nUnits
            count = count+1;
            xc_idx(:,count) = [i j];
        end
    end
    
    if 0 
        %%
        sps = 20000;
        xr = lags*binWidth/sps;
        mysort.plot.figure([1400 800])
        ah = subplot(1,2,1);
        plot(xr, xc(:,xc_idx(1,:) == xc_idx(2,:)))
        title('Autocorrelations');
        xlabel('time [sec]')
        
        ah(2) = subplot(1,2,2);
        plot(xr, xc(:,xc_idx(1,:) ~= xc_idx(2,:)))        
        title('Cross-Correlations');
        xlabel('time [sec]')
        set(ah, 'ylim', [0 1], 'xlim', xr([1 end]));
    end