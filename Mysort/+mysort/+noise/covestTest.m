%function covestTest()

    %% GENERATE DATA
    nC = 4;
    % Some noise
    X = 20*randn(nC,100000);
    % Heavily correlate first and last channel with a timeshift
    X(nC,:) = [X(1,3:end) [0 0]];
    % Correlate second and third with a timeshift
    X(2,:) = .5*X(2,:) + [X(3,2:end) [0]];
    
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
    
    % Decide for electrode positions (in this case 1D)
    el_positions = 1:size(X,1);
    % Build data handle for that configuration
    MEA = mysort.datasource.MultiElectrode(X, el_positions, 'srate', 32000);
    
    %% NOISE Estimation
    % Build noise estimator, dont care about template in the data
    C1 = mysort.noise.Covest(MEA, 'maxLag', Tf-1, 'forceMethod', 'matmul');
    C2 = mysort.noise.Covest(MEA, 'maxLag', Tf-1, 'forceMethod', 'xcorr');
    
    % Make a figure, showing the noise representation
    ccol = C1.buildCColumn(Tf+2);
    figure; imagesc(ccol);
    
    
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
    legend([h1 h2 h3], {'matched matmul', 'matched xcorr', 'template'});
   