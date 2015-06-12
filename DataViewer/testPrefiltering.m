X = DH.getData('trialIDs', 3902, 'tetrodeIDs', 36);

Y = util.prefilter(X);
chan = 1;

figure
plot(X(chan,1:200000))
hold on
plot(Y(chan,1:200000), 'r')

srate = 32000;

figure
mysort.util.frequencyAnalysis(X(chan,:), srate, 'onesided','.-');
hold on
mysort.util.frequencyAnalysis(Y(chan,:), srate, 'onesided','r.-');