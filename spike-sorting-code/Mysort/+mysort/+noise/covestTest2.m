function covestTest2()    
    channels1 = [1 10 20];% 30 40 50 60];
    Ls = 1000;
    maxLags = 10;
    [T1a T1b] = loop(maxLags,Ls,channels1);
    Ls = 10000;
    [T1a2 T1b2] = loop(maxLags,Ls,channels1);
    
    channels = 20;
    Ls2 = [1000 10000];% 50000 100000];
    maxLags = 10;
    [T2a T2b] = loop(maxLags,Ls2,channels);

    channels = 20;
    Ls = 10000;
    maxLags3 = [1 5 10 20 25 30 35 40 50];% 35 50 75 100];
    [T3a T3b T3c] = loop(maxLags3,Ls,channels);
    Ls = 20000;
    [T3a2 T3b2 T3c2] = loop(maxLags3,Ls,channels);
    

    mysort.plot.figure('w', 800, 'h', 900);
    ah = mysort.plot.subplot2([3 1], 'spacerY', 70, 'marginRight', 25);
    set(ah,'NextPlot', 'add');

    
    plot(ah(1), channels1, squeeze(T1a), '.-', 'linewidth', 2, 'markersize', 15);
    plot(ah(1), channels1, squeeze(T1a2), 'r.-', 'linewidth', 2, 'markersize', 15);
    plot(ah(1), channels1, squeeze(T1b), '.--', 'linewidth', 2, 'markersize', 15);
    plot(ah(1), channels1, squeeze(T1b2), 'r.--', 'linewidth', 2, 'markersize', 15);
    xlabel(ah(1), '# channels', 'fontsize', 14);
    legend(ah(1),'L=1000, lag=10', 'L=10000, lag=10','location', 'northwest');
    
    plot(ah(2), Ls2/1000, squeeze(T2a), '.-', 'linewidth', 2, 'markersize', 15);
    plot(ah(2), Ls2/1000, squeeze(T2b), '.--', 'linewidth', 2, 'markersize', 15);
    xlabel(ah(2), '# kSamples', 'fontsize', 14);
    legend(ah(2),'nC=20, lag=10', 'location', 'southeast');
    
    plot(ah(3), maxLags3, squeeze(T3a), '.-', 'linewidth', 2, 'markersize', 15);
    plot(ah(3), maxLags3, squeeze(T3b), '.--', 'linewidth', 2, 'markersize', 15);
    plot(ah(3), maxLags3, squeeze(T3c), '.--g', 'linewidth', 2, 'markersize', 15);
    xlabel(ah(3), 'maxlag', 'fontsize', 14);
    legend(ah(3),'nC=20, L=10000, xcorr', 'matmul', 'both', 'location', 'southeast');
    plot(ah(3), maxLags3, squeeze(T3a2), '.-k', 'linewidth', 2, 'markersize', 15);
    plot(ah(3), maxLags3, squeeze(T3b2), '.--k', 'linewidth', 2, 'markersize', 15);
    plot(ah(3), maxLags3, squeeze(T3c2), '.--c', 'linewidth', 2, 'markersize', 15);
    
    set(ah, 'fontsize', 14);
    ylabel('time [s]', 'fontsize', 14);
    %mysort.plot.savefig(gcf, 'runtime_analysis_new');

    function [Ta Tb Tc] = loop(maxLags, Ls, channels)
        Ta = zeros(length(maxLags), length(Ls), length(channels));
        Tb = Ta;
        Tc = Ta;
        for lag_idx=1:length(maxLags)
            lag = maxLags(lag_idx);
            for len_idx = 1:length(Ls)
                L = Ls(len_idx);
                for nC_idx = 1:length(channels)
                    nC = channels(nC_idx);

                    X = randn(nC, L);
                    el_positions = 1:nC;
                    MEA = mysort.ds.Matrix(X, 32000, el_positions);

                    %% NOISE Estimation
                    % Build noise estimators
                    tic;
                    C = mysort.noise.Covest(MEA, 'maxLag', lag, 'forceMethod', 'xcorr');
                    Ta(lag_idx, len_idx, nC_idx) = toc;
                    tic;
                    xc1 = C.xcovs;
                    C = mysort.noise.Covest(MEA, 'maxLag', lag, 'forceMethod', 'matmul');
                    Tb(lag_idx, len_idx, nC_idx) = toc;
                    tic;
                    xc2 = C.xcovs;
                    C = mysort.noise.Covest(MEA, 'maxLag', lag);
                    Tc(lag_idx, len_idx, nC_idx) = toc;                            
                    xc3 = C.xcovs;
                    tol = .0000001;
                    for i=1:size(xc1)
                        for k=1:size(xc1)
                            assert(all(xc1{i,k} <= xc2{i,k}+tol) && ...
                                   all(xc1{i,k} >= xc2{i,k}-tol), 'mismatch 1-2!');
                            assert(all(xc1{i,k} <= xc3{i,k}+tol) && ...
                                   all(xc1{i,k} >= xc3{i,k}-tol), 'mismatch 1-3!');
                        end
                    end
                end
            end
        end
    end
end