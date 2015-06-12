
nDimsRange = [1 2 5 10 50];
nR = length(nDimsRange);
nC = 1;
k = 3;
N = repmat(100, 1, k);
mysort.plot.figure('w', 800, 'h', 800);
D = []; ah = [];
for nDimsIdx = 1:length(nDimsRange);
    nDims = nDimsRange(nDimsIdx);
    mu = [0  zeros(1, nDims-1)
          5  zeros(1, nDims-1)
          10 zeros(1, nDims-1)];   
    
    [X P] = mysort.clustering.generateTestData('dim', nDims, 'k', k, 'mu', mu, 'N', N);
    for i=1:N(1)
        D(:,i) = sum( (X - repmat(X(i,:), size(X,1), 1)).^2, 2);
    end
    myGroup = D(1:N(1),:);
    others  = D(N(1)+1:end,:);
    ah(nDimsIdx) = subplot(nR, nC, (nDimsIdx-1)*1+1);
    hist(myGroup(:), 100)
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', 'g', 'EdgeColor','g')    
    hold on
%     subplot(nR, nC, (nDimsIdx-1)*1+1);
    hist(others(:), 100)
end
set(ah, 'xlim', [0 300])
