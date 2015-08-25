% function xcorrTest()
    st1 =    [10 20 30 100 110 120 130 200 300 300 301 302 400];
    st2 = 10+[10 20 30 100 110 120 130 200 300 300 301 302 400];
    L = max([st1 st2]);
    srate = 1000;
    IDs = [1 2];
    mysort.plot.xcorr({st1 st2}, 'srate', srate, 'IDs', IDs, 'T', L);
   

    rate = 30;
    duration = 800;
    st3 = poisson_spike_train(rate*2, duration);
    st4 = poisson_spike_train(rate, duration);
    st5 = sort([ st3(1:4:end)+.01 st4(1:10:end)-.01]); % poisson_spike_train(rate/10, duration)
    srate = 1000;
    st3 = round(st3*srate);
    st4 = round(st4*srate);
    st5 = round(st5*srate);
    L = duration*srate;
    mysort.plot.xcorr({st3 st4 st5}, 'srate', srate, 'IDs', [3 4 5], 'T', L);
    
    
%     assertEqual(idx, gtidx);