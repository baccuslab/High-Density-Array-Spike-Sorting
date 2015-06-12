maxSamples = 10000;
nEpochs = 50000;
nSamples = 1000000;
epochLength = 100;
maxLag = 4;
nC = 10;

cp = [1 1
      1 2
      1 3
      1 4
      2 1
      2 2
      2 3 
      2 4
      2 5
      2 6
      3 4
      3 6
      3 8
      4 4
      1 1
      2 2];
      
eStarts = floor(linspace(1, ceil(nSamples/2), nEpochs))';
epochs = [eStarts eStarts+repmat(epochLength-1, nEpochs, 1)]; 

X = randn(nSamples, nC);
X(2:end, 3) = X(2:end, 3) + X(1:end-1, 2);
X(:,4) = X(:,4)*2;
tic
xc = mysort.noise.computeXCorrs(X, cp, maxLag, epochs, maxSamples);
toc

