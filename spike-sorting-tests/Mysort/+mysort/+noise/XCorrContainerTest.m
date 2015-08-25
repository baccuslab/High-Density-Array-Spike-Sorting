%% INIT
nSamples = 100000;
nC = 10;
maxLag = 4;
spikeTrains = cell(1, nC);
for i=1:nC
    spikeTrains{i} = round(100 + (nSamples-200)*rand(10,1));
end
disp('generating data');
X = randn(nSamples, nC);

disp('correlating channel 2 and 3 with lag -1');
X(2:end, 3) = X(2:end, 3) + X(1:end-1, 2);

disp('doubling std of channel 4');
X(:,4) = X(:,4)*2;

disp('correlating channel 5 and 6 with lag 1');
X(2:end, 5) = X(2:end, 5) + X(1:end-1, 6);

%% Build XCorr container
XCC = mysort.noise.XCorrContainer(X, maxLag, 'spikeTrains', spikeTrains);

%% Invert example vector with levinson durbin
Tf = (maxLag+1);
inversion_channel_idx = [2 3 5 6];
nCinv = length(inversion_channel_idx);
x = ones(1,Tf*nCinv);
y = XCC.invMul(x, inversion_channel_idx);

figure; 
plot(x);
hold on
plot(y, 'g');

%% Build the different covariance represenations and plot them for various
%% channel sets
idxset = {1:4, 2:5, 3:6, 7:10, 1:10};
for i=1:length(idxset)
    idx = idxset{i};
    disp('Computation 1')
    tic
    xc = XCC.getXCorr4Channels(idx);
    toc
    disp('Computation 2')
    tic
    xc = XCC.getXCorr4Channels(idx);
    toc    
    ccol = XCC.getCCol4Channels(idx);
    Cte = XCC.getCte4Channels(idx);
    Cce = XCC.getCce4Channels(idx);
    
    

    figure;
    subplot(2,3,[1 4])
    imagesc(ccol)
    title('CCol for inversion')

    subplot(2,3,[2 3])
    imagesc(xc)
    title('internal xc functions')

    subplot(2,3,5)
    imagesc(Cte)
    title('C time embedded')

    subplot(2,3,6)
    imagesc(Cce)
    title('C channel embedded')
    mysort.plot.figureTitle(num2str(idx));
end