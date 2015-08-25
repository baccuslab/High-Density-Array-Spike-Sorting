nC = 10;
nS = 100;
X = randn(nS, nC)+10;

f1 = mysort.mea.filter_design(350, 8000, 20000, 10);
f2 = f1.copy();
mysort.mea.filter_init(f1, X);

Y1 = filter(f1, X);
Y2 = filter(f2, X);

spacer = mysort.plot.mc(X');
hold on
mysort.plot.mc(Y1', 'figure', 0, 'spacer', spacer, 'color', {'r'});
mysort.plot.mc(Y2', 'figure', 0, 'spacer', spacer, 'color', {'g'});
norm(Y1-X)
norm(Y2-X)

fvtool(f2)