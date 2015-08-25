X = randn(1000,10);
X(2:end-1,2) = X(2:end-1,2) + X(1:end-2,1) + X(3:end, 3);

X(10:end, 7) = X(1:end-9, 4);
xc = xcorr(X, 10, 'none');
Cte = mysort.noise.xcorr2Cte(xc);

figure;
subplot(2,1,1);
imagesc(xc);

subplot(2,1,2);
imagesc(Cte);