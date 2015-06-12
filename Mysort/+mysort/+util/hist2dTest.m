data = randn(12000,2);

[nc edgesx edgesy] = mysort.util.hist2d(data, 50, 50);

figure;
subplot(1,2,1)
plot(data(:,1), data(:,2), '.');
axis tight
subplot(1,2,2);
imagesc(nc)
axis normal xy
colorbar