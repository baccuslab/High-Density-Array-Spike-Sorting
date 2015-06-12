%testDistCluters

clear
profile on

nPtsPerClust = 2500;
nClust  = 3;
totalNumPts = nPtsPerClust*nClust;
m(:,1) = [1 1]';
m(:,2) = [-1 -1]';
m(:,3) = [1 -1]';
var = .6;
nDims = 2;
bandwidth = 1.3;
clustMed = [];
%clustCent;


x = var*randn(2,nPtsPerClust*nClust);
%*** build the point set
for i = 1:nClust
    x(:,1+(i-1)*nPtsPerClust:(i)*nPtsPerClust)       = x(:,1+(i-1)*nPtsPerClust:(i)*nPtsPerClust) + repmat(m(:,i),1,nPtsPerClust);   
end

% add some noise
x(:, end+1) = [-6 -1];
x(:, end+1) = [-6.4 -1.1];

x(:, end+1) = [-4   2.8-.8];
x(:, end+1) = [-4.4 2.7-.8];
x(:, end+1) = [-4.1 2.9-.8];
x(:, end+1) = [-4   3.4];
x(:, end+1) = [-4.2 3.3];

x(:, end+1) = [4 .9];
x(:, end+1) = [4 1.1];

%%
nMinPointsInCluster = 4;
maxNBandwidthIncreases = 1;
bandwidthIncreaseFactor = 1.25;
maxBWFact = 6;
cluFuns = {'Std MS', @(x) MeanShiftCluster(x, bandwidth)
           'InceaseBW', @(x) MeanShiftClusterIncreaseBW(x, bandwidth, 0 , nMinPointsInCluster, maxNBandwidthIncreases, bandwidthIncreaseFactor)
           'speedMS', @(x) MeanShiftClusterSpeed(x, bandwidth, maxBWFact)
           'BW+Speed', @(x) MeanShiftClusterSpeedBWInc(x, bandwidth, maxBWFact, nMinPointsInCluster, maxNBandwidthIncreases, bandwidthIncreaseFactor)
           };
nP = size(cluFuns,1);
mysort.plot.figure('w', 1000, 'h', 600);
for i=1:nP
    tic
    [clustCent,point2cluster,clustMembsCell] = cluFuns{i,2}(x);
    toc

    numClust = length(clustMembsCell);

    subplot(1, nP, i);
    hold on
    cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
    for k = 1:min(numClust,length(cVec))
        myMembers = clustMembsCell{k};
        myClustCen = clustCent(:,k);
        plot(x(1,myMembers),x(2,myMembers),[cVec(k) '.'])
        plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',10)
    end
    title([cluFuns{i,1} ' numClust:' int2str(numClust)])
end