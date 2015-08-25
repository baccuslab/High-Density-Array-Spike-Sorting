% function shiftMCRowsTest()    
    nC = 2;
    nS = 2;
    X = [ repmat([1:5 5:-2:1], nS, nC)];
    Y = mysort.util.shiftMCRows(X,[-1 2], nC);
    Z = mysort.util.shiftMCRows(Y,[-1 2], nC, 1);
    
    figure;
    mysort.plot.spikes(X, 'nC', nC);
    ax = gca; hold on
    mysort.plot.spikes(Y, 'nC', nC, 'ax', ax);
    
    