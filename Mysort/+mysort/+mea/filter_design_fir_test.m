
P.hpf = 300;
P.lpf = 7000;
sr = 20000;
P.fir_filterOrder = 110;
b  = mysort.mea.filter_design_fir(P.hpf, P.lpf, sr, P.fir_filterOrder);


Hd = dfilt.df2(b,1);
Hd.todf2tsos



X = randn(1000, 13);
X(100:200,1:2) = X(100:200,1:2) + 10;

Y1 = conv2(X, b(:), 'same');
Y2 = filtfilt(b(:),1,X);

figure
plot(Y1)
hold on
plot(Y2, ':')

