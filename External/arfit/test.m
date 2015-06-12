nC = 2;
C = eye(nC);
A = zeros(nC, 4);
w = zeros(1, nC);

X = arsim(w, A, C, 50000); X = X';
% mysort.plot.mc(X)

maxLag = 5;

Cest = mysort.noise.Covest(X, 'maxLag', maxLag);
Cest.getNoiseCovarianceMatrixTimeEmbed(2,[1:nC])

X2 = X; X2(2,:) = 2*X2(2,:);
Cest = mysort.noise.Covest(X2, 'maxLag', maxLag);
Cest.getNoiseCovarianceMatrixTimeEmbed(2,[1:nC])

c11 = [];
pred = [];
ar = 0:.1:1;

for i = 1:length(ar)
    a = ar(i);
    X3 = X; X3(1,2:end) = X(1,2:end)*(1-a) + X(2,1:end-1)*a;
    Cest = mysort.noise.Covest(X3, 'maxLag', maxLag);
    c = Cest.getNoiseCovarianceMatrixTimeEmbed(2,[1:nC]);
    c11(i) = c(1,1);
    a_ = [1-a a];
    pred(i) = a_*eye(2)*a_';
end
c
figure; plot(ar, c11, 'b.'); hold on; plot(ar, c11, 'g');

