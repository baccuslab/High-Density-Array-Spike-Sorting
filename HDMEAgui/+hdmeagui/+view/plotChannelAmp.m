function plotChannelAmp(H)
    ax = H.ChannelAmpAxes;
    cla(ax);
    axbg = H.ChannelAmpBGAxes;
    cla(axbg);

    G = H.GUI;
    
    linkaxes([ax axbg], 'xy');
    set(axbg, 'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'box', 'off')
    thr = 0.1;
    amps = G.WF.getFeatures('MinMax');
    amps = -amps(:,1);
    idx = amps(:,1) > thr;
    if ~any(idx)
        return
    end
    chans = G.WF.eventChans(idx);
    % plot brushable plots and link data
    h = plot(ax, chans, amps(idx), ...
        'kx', 'HitTest', 'off', 'markersize',4);
    maxi = max(amps(idx));
    set(ax, 'xlim', [0 G.WF.nC+1]);
    set(ax, 'ylim', [0 maxi*1.01]);
    xlabel(ax, 'channel index');
    ylabel(ax, 'amplitude [std noise]');
    set(ax, 'color', 'none');
    
    