% function covestTest3()
    % test also if the noise epoch setting works
    %% GENERATE DATA
    nC = 4;
    % Some noise
    L = 100000;
    X = 20*randn(nC,L);
        
    % Heavily correlate (make identical) first and last channel with a timeshift of 2
    X(nC,:) = [X(1,1:end) ]; %
    % Correlate second and third with a timeshift of +1
    X(2,:) = .5*X(2,:) + [    X(3,2:end) [0]];
    % Correlate first and fourth with a timeshift of -1
    %X(1,:) = .5*X(1,:) + [[0] X(4,1:end-1)  ];    
    %
%     figure;
%     plot(xcov(X(1,:), X(4,:)))
%     hold on
%     plot(xcov(X(4,:), X(1,:)), 'g')    
    
    X = X+.1*randn(nC,L);
    % Choose a "template"
    T = [-1 .5 -1 1
          0 0  0 0
          1 0  1 0
          -1 -2 -1 0];
    Tf = size(T,2);
    
    % Get channel representations
    t  = mysort.wf.m2v(T);
    t_ = mysort.util.embedTime2embedChan(t, nC);
    
    % Add template to data
    X(:,100+(1:Tf)) = X(:,100+(1:Tf))+T;
    X(:,400+(1:Tf)) = X(:,400+(1:Tf))+T;
    X(:,600+(1:Tf)) = X(:,600+(1:Tf))+T;
    X = X';
    noise_epochs = [620 L];
%     100
%                     120 400
%                     420 600
%                     620 L];
    % Decide for electrode positions (in this case 1D)
    el_positions = (1:size(X,2))';
    % Build data handle for that configuration
    MEA = mysort.ds.Matrix(X, 32000, 'test', el_positions);
    
    %% NOISE Estimation
    % Build noise estimator, dont care about template in the data
    C1 = mysort.noise.Covest2(MEA, 'maxLag', 15, ...
        'noiseEpochs', noise_epochs, ...
        'maxDist', 1000);
    C2 = mysort.noise.Covest2(MEA, 'maxLag', 15, ...
        'noiseEpochs', [], ...
        'maxDist', 1000);
    
    CC1 = mysort.noise.ccol2Cte(C1.CCol);
    CC2 = mysort.noise.ccol2Cte(C2.CCol);
    % Make a figure, showing the noise representation
    figure;
    subplot(1,2,1)
    imagesc(C1.CCol);
    subplot(1,2,2)
    imagesc(C2.CCol);   
 
    
    %% BUILD FILTER
    % use the usual equation f = inv(C)*t / sqrt( t' * inv(C) * t)
    f1_ = C1.invMul(t_); %f1_ = f1_/sqrt(t_*f1_');
    f2_ = C2.invMul(t_); %f2_ = f2_/sqrt(t_*f2_');
    % get filter representations
    f1  = mysort.util.embedChan2embedTime(f1_,nC);
    f2  = mysort.util.embedChan2embedTime(f2_,nC);
    F1  = mysort.wf.v2m(f1, nC)
    F2  = mysort.wf.v2m(f2, nC)
    
    %% FILTER DATA
    % with the matched filter
    y1 = mysort.util.mcfilt(X, F1);
    y1(1:Tf) = 0;
    y2 = mysort.util.mcfilt(X, F2);
    y2(1:Tf) = 0;    
    % with the template
    z = mysort.util.mcfilt(X, T);
    
    %% PLOT DATA
    spacer = mysort.plot.mc(X(:,1:1000), 'color', {'k'});
    hold on
    h1 = plot(y1(1:1000), 'g');
    h2 = plot(y2(1:1000), 'g:');
    h3 = plot(z(1:1000), 'r');
    legend([h1 h2 h3], {'noise epoch', 'full', 'template'});
   