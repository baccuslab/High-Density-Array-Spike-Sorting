function block_levinson_test(C, xi, nC)
if ~exist('nC', 'var')
    nC = 1;
end

f = inv(C)*xi;
g = C\xi;

d = size(C,1)/nC;

l = matlabfilecentral.block_levinson(xi, C(:,1:d));
f = f/norm(f);
g = g/norm(g);
l = l/norm(l);

norm(l-f)

F = mysort.util.v2m(f',nC);
G = mysort.util.v2m(g',nC);
L = mysort.util.v2m(l',nC);
X = mysort.util.v2m(xi',nC);


figure;
subplot(2,1,1)
plot(f, 'marker', '.'); hold on; plot(l,'r'); plot(g,'g');
title('filter coefficients');

subplot(2,1,2)
plot(mysort.util.mcfilt(X,F), 'marker', '.');
hold on
% plot(mysort.util.mcfilt(X,fliplr(F)), 'b:');

plot(mysort.util.mcfilt(X,L), 'r');
% plot(mysort.util.mcfilt(X,fliplr(L)), 'r:');
plot(mysort.util.mcfilt(X,G), 'g');
title('convolution functions');
legend('normal inversion', 'levinson inversion', 'C\\xi');