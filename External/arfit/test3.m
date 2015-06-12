% Define AR model
nC = 2;
A = [.7 .2 .0 .0
     .2 .7 .1 .0];
%B = [A; eye(nC) zeros(nC)];
c = [3  0
     0  1];
w = zeros(1, nC);

% Simulate data
X = arsim(w, A, c, 500000); X = X';
% mysort.plot.mc(X(:, 1:10000));

% Fit AR model
pmin = 1;
pmax = 3;
[w_, A_, c_, sbc_, fpe_, th_] = arfit(X', pmin, pmax);

% Calculate Covariance matrix
maxLag = 5;
Cest = mysort.noise.Covest(X, 'maxLag', maxLag);
C = Cest.getNoiseCovarianceMatrixTimeEmbed(5,[1:nC]);

% Derive AR model from Covariance matrix
order = 3;
[A__, gamma_y, c__] = mysort.noise.AR_from_C(C, nC, 2);
A
A_
A__
c
c_
c__


