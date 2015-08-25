
function [C H] = plot2(h, lims, x, noise_idx, spike_idx, tau, i1, i2, blegend)
axes(h);
if blegend
    plot(1000,1000, '.r', 'markersize', 12);
    hold on
    plot(1000,1000, '.k', 'markersize', 12);
    
    H = mysort.noise.plot(h, lims, x, spike_idx, tau, i1, i2, 'r.');
    C = mysort.noise.plot(h, lims, x, noise_idx, tau, i1, i2);
    
    legend('spikes', 'noise', 'Location', 'SouthEast');
    legend('boxoff')
else
    H = mysort.noise.plot(h, lims, x, spike_idx, tau, i1, i2, 'r.');
    hold on
    C = mysort.noise.plot(h, lims, x, noise_idx, tau, i1, i2, 'k.');
end

