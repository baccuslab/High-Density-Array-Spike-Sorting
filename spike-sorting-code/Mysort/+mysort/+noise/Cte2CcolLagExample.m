
nC = 10;
X = randn(1000,nC);
X(2:end-1,2) = X(2:end-1,2) + X(1:end-2,1) + X(3:end, 3);

X(10:end, 7) = X(1:end-9, 4);
xc = xcorr(X, 10, 'none');
Cte = mysort.noise.xcorr2Cte(xc);

CcolLag0 = mysort.noise.Cte2CcolLag(Cte, nC, 0);
CcolLag1 = mysort.noise.Cte2CcolLag(Cte, nC, 1);
CcolLag9 = mysort.noise.Cte2CcolLag(Cte, nC, 9);

figure;
subplot(1,5,1);
imagesc(xc);
title('Matlab xcorr');

subplot(1,5,2);
imagesc(Cte);
title('Time embedding');

subplot(1,5,3);
imagesc(CcolLag0);
title('Lag 0');

subplot(1,5,4);
imagesc(CcolLag1);
title('Lag 1');

subplot(1,5,5);
imagesc(CcolLag9);
title('Lag 9');