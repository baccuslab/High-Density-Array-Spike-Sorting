
mu      = [.5 .5 0]*sqrt(2);
mu(2,:) = [ 0  0 1];

[nMu dims] = size(mu);
nS = 10000;

X = zeros(nS, dims);
for i=1:nMu
    X = X + repmat(mu(i,:), nS, 1) .* repmat(randn(nS,1), 1, dims);
end

noiseStd = 0.1;
N = noiseStd*randn(nS, dims);

Y=X+N;

figure;
plot3(X(:,1), X(:,2), X(:,3), '.');

C = cov(X);

[D V] = eig(C)

cvar = max(diag(C));
Q = C/cvar;
% RML = eye(dims) - X*inv(X*inv(Q)*X')*X'*inv(Q);