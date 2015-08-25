figure;
a = axes();

P = struct();
P.axesHandles = a;
P.plotNormFit = 1;
P.edges = [];
P.oversampleFit = 10;
P.normFitPara = {'color', 'r', 'linewidth', 2, 'marker', 'o'};


x = randn(1000,1);

mysort.plot.normalityHist(x, P); 