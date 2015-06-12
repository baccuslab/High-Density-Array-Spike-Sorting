function alignWaveformsUpsampleMaxOrMinTest()
warning('This function is depricated! Use mysort.wf.* instead!');
    nS = 100;
    nC = 2;

    X = repmat([5:-.5:0 1:10 10:-1:1 0:.5:5], nS, nC)-5;
    X = X + randn(size(X));
    shifts = [-4:2:4];
    stau = repmat(shifts', nS/length(shifts),1);
    X = mysort.util.shiftMCRows(X,stau,nC,1);
    debug = 1;
    mysort.plot.spikes(X,'nC',nC);
    title('Not aligned Spikes'); 
%     [tau Y] = mysort.util.alignWaveformsOnMax(X, 2, 'debug',1);
% 
%     fprintf('Shift Error after Iteration %d: %d\n', i, sum(abs(-tau-stau)));
%     mysort.plot.spikes(Y,'nC',nC);
%     title('Aligned Spikes');
    
    [tau Y] = mysort.util.alignWaveformsUpsampleMaxOrMin(X, 2, 'maxIdx', 20);
    fprintf('Shift Error after Iteration %d: %.4f\n', i, sum(abs(-tau-stau)));
    mysort.plot.spikes(Y,'nC',nC);
    title('Aligned Spikes');    
%     [tau Y] = mysort.util.alignWaveformsOnMax(X, 2, 'debug',1, 'maxIdx', 30);
%     fprintf('Shift Error after Iteration %d: %d\n', i, sum(abs(-tau-stau)));
%     mysort.plot.spikes(Y,'nC',nC);
%     title('Aligned Spikes');        

