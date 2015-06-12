close all;
clear all
nC = 10;
Tf = 35;
nE = 1000;
X = randn(nE, nC*Tf);

WF = mysort.wf.WfMatrix(X, 1:nE, ones(1, nE), 1:nE, 10, Tf);

wfp = mysort.plot.WaveformManagerPlot(WF, 'plotMedianColor', [.9 .7 .0]);