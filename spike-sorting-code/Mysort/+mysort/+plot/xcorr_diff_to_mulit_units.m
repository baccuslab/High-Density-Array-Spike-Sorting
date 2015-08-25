
function [xc_unnormalized bin_centers] = xcorr_diff_to_mulit_units(spike_trains, srate, maxlag, trialLength)

    mua_others = mysort.spiketrain.leave_one_out_multi_units(spike_trains);

    xc = []; xc_unnormalized=[];
    for unit = 1:size(spike_trains,2)
        [xc(unit,:), E_xc, var_xc, conf_95_xc, bin_centers, edges, xc_unnormalized(unit,:)] = ...
            mysort.spiketrain.xcorr(spike_trains(:,unit), mua_others(:,unit), ...
            'T', trialLength, 'maxLag', maxlag);
    end
    mysort.plot.figure('name', 'XCorr Diff to Mulitunit');
    ah = mysort.plot.subplot(size(spike_trains,2));
    for unit = 1:size(spike_trains,2)
        plot(ah(unit), bin_centers/srate*1000, xc_unnormalized(unit,:), 'k');
    end
