L = 100000;
X = randn(L,1);
X = filter(ones(1,5)/5,1,X);
X = X/var(X);
t = 2*[0 1 2 1 -1 -2 -1 0 .5 1 0]';
spikeTimes = (100:400:L)';
nS = length(spikeTimes);
gdf = [ones(nS, 1) spikeTimes];
[x, y] = mysort.wf.templateSpikeTrainData(t, gdf, 0);
idx = x(~isnan(x));
X(idx) = X(idx) + y(~isnan(x))';

% xc = xcorr(X, maxLag, 'none')';
% figure, plot(xc)
% ac = xc(maxLag+1:end);
% figure, plot(ac)
% H = toeplitz(ac);
% figure, imagesc(H)

Tf = 21;
Tf2 = floor(Tf/2);
maxLag = Tf-1;
Hest = mysort.noise.Covest2(X, 'maxLag', maxLag);
H = mysort.noise.ccol2Cte(Hest.CCol);
iH = inv(H);
% figure, imagesc(H)
% figure, imagesc(iH)

nIter = 20;
L = size(X,1);
f = randn(Tf, 1);
d = 2;
for i=1:nIter
    i
    y = conv2(X, f, 'same');
    assert(~any(isnan(y)), 'nan!')
    [pks, locs] = findpeaks(y, 'MINPEAKDISTANCE', Tf, 'NPEAKS', 200);
    invalid = locs-Tf2 <=0 | locs+Tf2 > L;
    pks(invalid) = [];
    locs(invalid) = [];
    
    N = 0;
    w = pks.^d;
    T = zeros(Tf, length(locs));
    for k=1:length(locs)
        T(:,k) = w(k)*X(locs(k):locs(k)+Tf-1);
    end
    xi =sum(T, 2)/sum(w);
    if 0
        %%
        figure
        plot(F)
        hold on
        plot(xi, 'k', 'linewidth', 3)
%         plot(t, '--g', 'linewidth', 2)
    end
    
    f = iH * xi;
    n = (f'*H*f);
    assert(n>0, '0!')
    f = f/ n;    
    assert(~any(isnan(f)), 'nan!')
end
%%
% k = randn(Tf,1);
% y = conv2(X, k, 'same');
% var(y)
% k'*H*k
% k = k/sqrt(k'*H*k);
% y = conv2(X, k, 'same');
% var(y)
% k'*H*k

%%
figure
subplot(4,1,1)
plot([X y], 'linewidth', 2)
title('Data and filter output');

subplot(4,1,2)
plot(t)
title('original template');

subplot(4,1,3)
plot(xi)
title('estimated template');

subplot(4,1,4)
plot(f)
title('estimated filter');