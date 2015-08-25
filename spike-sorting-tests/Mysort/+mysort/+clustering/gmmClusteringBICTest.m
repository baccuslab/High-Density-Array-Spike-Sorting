%function gmmClusteringBICTest()

nDims = 2;
[X P] = mysort.clustering.generateTestData('dim', nDims, 'k', 10);
% nDims = 2;
% k = 5;
% N = [150 117 227 137 134];
% [X P T] = mysort.clustering.generateTestData('dim', nDims, 'k', k, 'mu_spread', 5, 'N', N);
% 
% 
[c centers N obj R] = mysort.clustering.gmmClustering(X, ...
    'repeats', 50, ...
    'kmin', 1, 'kmax', 15, 'prewhitened', true, 'homoscedastic', true,...
    'fixedCov', ones(1,nDims)); %[]);


rangeX = [min(X(:,1)) max(X(:,1))];
rangeY = [min(X(:,2)) max(X(:,2))];

mysort.plot.figure('w', 700, 'h', 700);
scatter(X(:,1),X(:,2),10,'.')
hold on
pdffun  = @(x,y) pdf(obj,[x y]);
plotfun = @(x,y) bsxfun(@max, pdffun(x,y), 0.0001);
plotfun = @(x,y) log(pdffun(x,y)-.000001);
h = ezcontour(plotfun, rangeX, rangeY, 100);

plot(P.mu(:,1), P.mu(:,2), 'rx', 'markersize', 20, 'linewidth', 4);
plot(centers(:,1), centers(:,2), 'ro', 'markersize', 20, 'linewidth', 4);

fprintf('# Cluster: %d, min(AIC): %d, min(BIC): %d\n', ...
    P.k, R.nClusterAIC, R.nClusterBIC);
fprintf('Elem in Cluster: ');
fprintf('%5d ', P.N);
fprintf('\nElem clustered:  ');
fprintf('%5d ', N);
fprintf('\n');


