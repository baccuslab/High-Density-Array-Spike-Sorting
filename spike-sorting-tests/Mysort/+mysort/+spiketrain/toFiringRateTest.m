
fs = .01;
w  = .1; 
L = 1;
ST = [ .5 .55 .575 .58 .8 .85 .1:.001:.3];
[fr x] = mysort.spiketrain.toFiringRate(ST, w, fs, [0 L]);

figure
ah = subplot(2,1,1);
mysort.plot.spiketrain({ST})
ah(2) = subplot(2,1,2);
plot(x, fr)

linkaxes(ah, 'x');


%%
fs = .01;
w  = .1; 
L = 1;
ST = linspace(0, 1, 1000);
[fr x] = mysort.spiketrain.toFiringRate(ST, w, fs, [0 L]);

figure
ah = subplot(2,1,1);
mysort.plot.spiketrain({ST})
ah(2) = subplot(2,1,2);
plot(x, fr)

linkaxes(ah, 'x');