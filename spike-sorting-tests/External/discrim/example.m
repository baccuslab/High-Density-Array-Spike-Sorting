[X P] = mysort.clustering.generateTestData();

r = randn(3);
C = r'*r;     % random positive definite symmetric matrix
M = randn(4, 3)*2;    % random means
k = ceil(rand(400, 1)*4);    % random class indeces
X = randn(400, 3)*chol(C) + M(k, :);  
f = lda(X, k); disp(f)
[lambda ratio] = cvar(f)   % canonical variates
cov(f), plotcov(f)
plotcov(shrink(f, .5))
[c post] = classify(f, X);
confmat(k, c)
confmat(k, post)