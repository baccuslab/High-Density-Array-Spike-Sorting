% function alignWaveformsInterpolatedMeanTest()
warning('This function is depricated! Use mysort.wf.* instead!');
    nS = 100;
    nC = 2;
    maxShift = 9;
    X = repmat([zeros(1, maxShift) [5:-.5:0 1:10 11 10:-1:1 0:.5:5]-5 zeros(1, maxShift) ...
                zeros(1, maxShift) -[5:-.5:0 1:10 11 10:-1:1 0:.5:5]+5 zeros(1, maxShift)], nS, nC/2);
    mysort.plot.spikes(X,'nC',nC);
    title('raw traces');
    
    X = X + randn(size(X));
    mysort.plot.spikes(X,'nC',nC);
    title('raw traces with noise');
    
    shifts = [-maxShift:2:maxShift];
      
%     X = repmat([5:-.5:0 1:10 10:-1:1 0:.5:5], nS, nC)-5;
%     X = X + randn(size(X));
%     shifts = [-4:2:4];
    
    stau = repmat(shifts', nS/length(shifts),1);
    X = mysort.util.shiftMCRows(X,stau,nC,1);
    debug = 1;
    mysort.plot.spikes(X,'nC',nC);
    title('shifted traces');  
    
%     [tau Y] = mysort.util.alignWaveformsOnMax(X, 2, 'debug',1);
% 
%     fprintf('Shift Error after Iteration %d: %d\n', i, sum(abs(-tau-stau)));
%     mysort.plot.spikes(Y,'nC',nC);
%     title('Aligned Spikes');
    
    [tau Y] = mysort.util.alignWaveformsInterpolatedMean(X, 2);
    fprintf('Shift Error after Alignment %.4f: \n', sum(abs(-tau-stau)));
    mysort.plot.spikes(Y,'nC',nC);
    title('Aligned Traces Mean');    
   
    [tau_min Y] = mysort.util.alignWaveformsInterpolatedMin(X, 2);
    fprintf('Shift Error after Alignment %.4f: \n', sum(abs(-tau_min-stau)));
    mysort.plot.spikes(Y,'nC',nC);
    title('Aligned Traces Min');    

