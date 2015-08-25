% function xcorrTest2()
gdf = gdfbotmcleaned(ismember(gdfbotmcleaned(:,1),[1]),:);
mysort.plot.xcorr(gdf(1:10000,:), 'maxLag', 150)

%%
    rate = 30;
    duration = 800;
    srate = 32000;
    L = duration*srate;
   
    
    st3 = poisson_spike_train(rate*2, duration);
    st4 = poisson_spike_train(rate, duration);
    st5 = sort([ st3(1:4:end)+.01 st4(1:10:end)-.01]); % poisson_spike_train(rate/10, duration)
    st3 = round(st3*srate);
    st4 = round(st4*srate);
    st5 = round(st5*srate);
    
    
    mysort.plot.xcorr({st3 st4 st5}, 'srate', srate, 'IDs', [3 4 5], 'T', L);
    mysort.plot.xcorr({st3 st4 st5}, 'srate', srate, 'IDs', [3 4 5], 'T', L, 'binSize', 1, 'maxLag', 500);
    
    
%     assertEqual(idx, gtidx);

%%
    rate2 = 100;
    rate1 = 1;
    duration = 80;
    srate = 32000;
    L = duration*srate;    
    
    st3 = poisson_spike_train(rate1, duration);
    st4 = poisson_spike_train(rate2, duration);
    st5 = sort([ st3(1:4:end)+.01 st4(1:10:end)-.01]); % poisson_spike_train(rate/10, duration)
    st3 = round(st3*srate);
    st4 = round(st4*srate);
    st5 = round(st5*srate);
    
%     mysort.plot.xcorr({st3 st4 st5}, 'srate', srate, 'IDs', [3 4 5], 'T', L);
    mysort.plot.xcorr({st3 st4 st5}, 'srate', srate, 'IDs', [3 4 5], 'T', L, 'binSize', 1, 'maxLag', 500);
    
%%
    rate1 = 1;
    rate2 = 100;
    duration1 = 20;
    duration2 = 20;
    srate = 32000;
    L = (duration1+duration2)*srate;    
    
    st1 = poisson_spike_train(rate1, duration1);
    st2 = poisson_spike_train(rate2, duration2);
    
    st = round([st1 st2+duration1]*srate);
    mysort.plot.xcorr([ones(length(st),1) st'], 'srate', srate, 'IDs', [1], 'T', L, 'binSize', 1, 'maxLag', 500);