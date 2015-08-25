%% INIT "true" AR model
nC = 2;
A = [.5 .2 .1 0
     .2 .5 0  .1];
w = [0 0];
c = [1 .5
    .5 1 ];

%% Simulate data with the "true" model
X = arsim(w, A, c, 50000); X = X';
pmin = 1; pmax = 3;

%% Fit model to simulated data
[w_, A_, c_, sbc_, fpe_, th_] = arfit(X', pmin, pmax);

%% Compute Cov on simulated data
maxLag = 5;
Cest = mysort.noise.Covest(X, 'maxLag', maxLag);
C = Cest.getNoiseCovarianceMatrixTimeEmbed(4,[1:nC]);

%% Get AR model from cov
[A__, gamma_y, c__] = mysort.noise.AR_from_C(C, nC, 3);

%% Check for consitency
A
A_
A__
c 
c_
c__
