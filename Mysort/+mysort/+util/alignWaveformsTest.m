warning('This function is depricated! Use mysort.wf.* instead!');
close all
clear all
nS = 100;
nC = 2;

X = repmat([5:-.5:0 1:10 9:-1:1 0:.5:5], nS, nC)-5;
X = X + 1*randn(size(X));
shifts = [-4:2:4];
stau = repmat(shifts', nS/length(shifts),1);
X = mysort.util.shiftRows(X,stau,nC,1);
debug = 1;
nettoShift = zeros(size(stau));

tau = mysort.util.alignWaveforms(X, 2, 'debug',0, 'nIter',4);
X = mysort.util.shiftRows(X,tau,nC,1);

fprintf('Shift Error: %d\n',  sum(abs(-tau-stau)));
mysort.plot.spikes(X,'nC',nC);
title('Aligned Spikes');

