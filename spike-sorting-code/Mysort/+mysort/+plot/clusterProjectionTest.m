K = 3;
N = 200;
X = [];
D = 10;
M = [];
c = [];
for i=1:K
    M(i,:) = i*randn(1,D);
    X = [X; randn(N, D) + repmat(M(i,:),N,1)];
    c = [c; ones(N,1)*i];
end

mysort.plot.clusterProjection(X, c);