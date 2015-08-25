
nC = 10;
X = randn(1000,nC);
X(2:end-1,2) = X(2:end-1,2) + X(1:end-2,1) + X(3:end, 3);

X(10:end, 7) = X(1:end-9, 4);
xc = xcorr(X, 10, 'none');
Cte = mysort.noise.xcorr2Cte(xc);

Cce = mysort.noise.Cte2Cce(Cte, nC);

figure;
subplot(1,3,1);
imagesc(xc);
title('Matlab xcorr');

subplot(1,3,2);
imagesc(Cte);
title('Time embedding');

subplot(1,3,3);
imagesc(Cce);
title('Channel embedding');