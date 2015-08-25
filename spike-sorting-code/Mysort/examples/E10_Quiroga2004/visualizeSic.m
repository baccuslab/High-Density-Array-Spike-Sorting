if 1
    path_to_simulated_files = E10_quirogaDataPath();

    benchmarks = {'Easy1', 'Easy2', 'Difficult1', 'Difficult2'};
    b = 1; noisestr = '015';
    s1 = 114174; s2 = s1+10000; % this is a sorting error
    s1 = 152140; s2 = s1+600;
    s1 = 76140; s2 = s1+300;
    quirogaFile = sprintf('C_%s_noise%s.mat', benchmarks{b}, noisestr);

    Tf = 71; cutLeft = -10;
end
if 0
    GT = E10_preprocessing(path_to_simulated_files, quirogaFile, cutLeft, Tf);
end
%%
if 1
    T = GT.templates;

    C = GT.C + eye(size(GT.C))*GT.C(1,1)/8;

    F = T/C;
    nF = size(F,1);

    upsample = 3;
    spikePrior = .01;


    X = GT.X(:, s1:s2);
    gdf = GT.gdf(GT.gdf(:,2)>=s1 & GT.gdf(:,2)<=s2,:);
    gdf(:,2) = gdf(:,2) - s1;
    gdf(gdf(:,1)==2,2) = gdf(gdf(:,1)==2,2)+5;
    gdf(gdf(:,1)==3,2) = gdf(gdf(:,1)==3,2)+1;
    nC = size(X,1);
    noisePrior = 1-nF * spikePrior; 
    threshold = log(noisePrior);
               % Calculate confusion matrix
    CONF    = mysort.util.calculateXIvsF(T, F, nC, 0);
    CONF_up = mysort.util.resampleTensor(CONF, upsample, 1);

    Tf_up = Tf * upsample;
    blockLen = ceil(Tf_up/4);  
    priors = repmat(spikePrior, nF, 1); 
end

% mysort.plot.sorting(GT.X, GT.gdf, GT.templates);

% close all
Y = mysort.util.applyMcFilters(X, F);
[D Dshifts] = mysort.util.calculateDiscriminantFunctions(Y, CONF, priors);

[mPr c] = max(D,[],1);
epochs = mysort.epoch.fromBinaryVectorMinLen(...
                        mPr>threshold, round(1.2*Tf));
nE = size(epochs,1);      
R = struct;
for e=1:nE
    R(e).D = [];
    R(e).S = [];
    R(e).p = [];
    if upsample>1
        Yup = [];
        % UPSAMPLE THE FILTER OUTPUTS NOT THE DISCRIMINANT
        % FUNCTIONS !!! D will be not zero mean and you will get
        % resampling artifacts at the right end!
        Yup(1:nF,:) = mysort.util.resampleMC(...
            Y(1:nF,epochs(e,1):epochs(e,2)), upsample, 1);
    else
        Yup = Y(:,epochs(e,1):epochs(e,2));
    end
    Dup = mysort.util.calculateDiscriminantFunctions(Yup, CONF_up, priors);
    [maxD maxClasses] = max(Dup,[],1); 
    count = 1;
    R(e).D{1} = Dup;
    while any(maxD>threshold)
        [maxd t] = max(maxD,[],2);
        maxC = maxClasses(t);
        offset = t - (Tf_up - (upsample-1));
        subtractor = zeros(size(Dup));
        for f=1:nF
            sub = squeeze(CONF_up(:,f,maxC))';
            subtractor(f,:) = mysort.util.shiftSubtract(zeros(1,size(Dup,2)),...
                                sub, offset, false);

        end
        Dup = Dup + subtractor + log(priors(maxC));
        R(e).D{count+1} = Dup;
        R(e).S{count+1} = subtractor;
        R(e).p(count+1) = log(priors(maxC));
        [maxD maxClasses] = max(Dup,[],1);                   
        count = count +1;
    end
end
%%     


for e=1:nE
    f = mysort.plot.figure('w', 1000, 'h', 500);
    nP = 2 + (length(R(e).S)-1)*2;
    p = 1; ah = zeros(nP,1);
    epochs(e,1)
    % plot data
    x = X(:, epochs(e,1):epochs(e,2));
    gdf_e = gdf(gdf(:,2)>=epochs(e,1) & gdf(:,2)<=epochs(e,2),:);
    gdf_e(:,2) = gdf_e(:,2) - epochs(e,1);
    title('Data with ground truth', 'fontsize', 14);
    ylabel('Voltage [a.u.]', 'fontsize', 14); 
   
    ah(p) = subplot(nP,1,p);
    mysort.plot.sorting(x, gdf_e, T, 'axesHandle', ah(p));
    p=p+1;
    
    % plot D
    ah(p) = subplot(nP,1,p);
    plot([0 size(R(e).D{1},2)], threshold*[1 1], ':', 'color', [.4 .4 .4], 'linewidth', 2);
    hold on
    D = R(e).D{1};
    for k=1:3
        plot(D(k,:), 'linewidth', 2, 'color', mysort.plot.vectorColor(k));  
    end
    title('Discriminant functions (iteration 0)', 'fontsize', 14);
    p=p+1;
    
    for it = 2:length(R(e).S)
        % plot subtractor
        ah(p) = subplot(nP,1,p);
        hold on
        SUB = R(e).S{it} + R(e).p(it);
        for k=1:3
            plot(SUB(k,:), 'linewidth', 2, 'color', mysort.plot.vectorColor(k));  
        end    
        title(sprintf('Subtractor (iteration %d)', it-1), 'fontsize', 14);
        p=p+1;

        % plot D2
        ah(p) = subplot(nP,1,p);
        plot([0 size(R(e).D{1},2)], threshold*[1 1], ':', 'color', [.4 .4 .4], 'linewidth', 2);
        hold on
        D = R(e).D{it};
        for k=1:3
            plot(D(k,:), 'linewidth', 2, 'color', mysort.plot.vectorColor(k));  
        end
        title(sprintf('Discriminant functions (iteration %d)', it-1), 'fontsize', 14);
        p=p+1;
    end
    
    for p=1:nP
        axes(ah(p));
        set(ah(p), 'fontsize', 14);
        box off
        if p<nP
            set(ah(p), 'xtick', [], 'xticklabel', []);
        else
            xlabel('time [samples]', 'fontsize', 14);
        end
        axis tight
    end    
    linkaxes(ah([2 4]), 'xy');
    linkaxes(ah([2:4]), 'x');

end