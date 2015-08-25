nC = 100;

%%
L = 10^3;
X = randn(L,nC);
tic
xcovs = xcorr(X, 60);
toc
% 1.35 sec


%%
error('this produces already out of memory!');
L = 10^4;
X = randn(L,nC);
tic
xcovs = xcorr(X, 60);
toc