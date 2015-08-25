x = 0:.1:30;
w = x(2)-x(1);
dims = 10;
t = randn(1,10);

e = norm(t);
scalefactor = 5;
y1 = ncx2pdf(x, dims, e^2);
y2 = scalefactor^2*ncx2pdf(scalefactor^2*x, dims, e^2);
% y1 = chi2pdf(x, dims);
% y2 = chi2pdf(scalefactor^2*x, dims);

 
bS = 10000;
X = repmat(t, nS, 1) + randn(nS,dims);
Y1 = histc(sum(X.^2,2), x)/(w*nS);
Y2 = histc(sum((X/scalefactor).^2,2), x)/(w*nS);



figure;
subplot(2,1,1)
bar(x+w/2,Y1)
hold on
plot(x, y1, 'g', 'linewidth',2);

subplot(2,1,2)
bar(x+w/2, Y2)
hold on
plot(x, y2, 'g', 'linewidth',2);
set(gca, 'xlim', [0 3])
