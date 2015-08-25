function P = xcorrBinned(spiketrain_or_gdf, binWidth, varargin)
P.sps = 20000;
P.XCS = [];
P = mysort.util.parseInputs(P, varargin, 'error');

if isempty(P.XCS)
    P.XCS.binWidth = binWidth;
    [P.XCS.xc, P.XCS.edges, P.XCS.binnedST, P.XCS.lags, P.XCS.xc_idx] = mysort.spiketrain.xcorrBinned(spiketrain_or_gdf, binWidth, varargin{:});
end

xr = P.XCS.lags*P.XCS.binWidth/P.sps;
auto = P.XCS.xc(:,P.XCS.xc_idx(1,:) == P.XCS.xc_idx(2,:));
ccs  = P.XCS.xc(:,P.XCS.xc_idx(1,:) >  P.XCS.xc_idx(2,:));

P.fh = mysort.plot.figure([1400 800]);
P.ah = subplot(1,2,1);

auto((size(auto,1)+1)/2,:) = nan;
plot(xr, auto)
title('Autocorrelations (0-bin masked)');
xlabel('time [sec]')

P.ah(2) = subplot(1,2,2);
plot(xr, ccs)        
title('Cross-Correlations');
xlabel('time [sec]')
set(ah, 'ylim', [0 1], 'xlim', xr([1 end]));