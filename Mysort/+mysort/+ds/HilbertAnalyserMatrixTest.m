nSample = 10000;
nChannel = 10;
X = randn(nSample, nChannel);
samplesPerSecond = 20000;
hpf = 800;
lpf = 1200;
order = 6;
H = mysort.ds.HilbertAnalyserMatrix(X, samplesPerSecond, hpf, lpf, order);
F = mysort.ds.FrequencyAnalyserMatrix(X, samplesPerSecond, hpf, lpf, order);
mysort.plot.SliderDataAxes({X F H}, 'channelSpacers', [5 .5 0]);