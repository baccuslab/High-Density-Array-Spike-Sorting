




%%







%%
figure;
hist(smad);
% figure;





%%

figure;
hist(D, 500);

%%
figure;
h=mysort.plot.subplot2([5,1], 'max');
plot(Y(nC,:), 'k');
hold on
mysort.plot.epochs(NE{ch}, 0, 'c', 'linewidth', 5);
mysort.plot.epochs(NEN{ch}, 0, 'm', 'linewidth', 3);
for n=1:4%length(neighbors)
    axes(h(n+1));
    ni = neighbors(n);
    plot(Y(ni,:),'k')
    hold on
    mysort.plot.epochs(NE{ni}, 0, 'c', 'linewidth', 5);
    mysort.plot.epochs(NEN{ni}, 0, 'm', 'linewidth', 3);
end
linkaxes(h,'xy');

%%
figure;
max_xc = max(xcovs, [], 2);
plot(D,max_xc, '.');

figure;
hist3([D' max_xc])
%%
Tf = (size(xcovs,2)+1)/2;
lagrange = -Tf+1:Tf-1;
idx = find(D<200&D>150, 1);

figure; plot(lagrange, xcovs(idx,:))

%%
chan1 = IDX(idx,1);
chan2 = IDX(idx,2);

figure;
h1=subplot(2,1,1);
plot(Y(chan1,:),'k');
hold on
comnoise = mysort.epoch.intersect(NE{chan1},NE{chan2});
mysort.plot.epochs(comnoise, 0, 'r', 'linewidth',4);
h2=subplot(2,1,2);
plot(Y(chan2,:),'k');
hold on
mysort.plot.epochs(comnoise, 0, 'r', 'linewidth',4);
linkaxes([h1 h2],'xy');


%%
eidx = 10;
neighbors = mysort.mea.nearestElectrodes(channel_x, channel_y, eidx, 6);
mysort.plot.figure;
for i=1:length(channel_x)
    plot(channel_x(i), channel_y(i), '.k');
    hold on
end
plot(channel_x(eidx), channel_y(eidx), 'or');
hold on
for i=1:length(neighbors)
    eidx = neighbors(i);
    plot(channel_x(eidx), channel_y(eidx), 'ob');
end

%%
eidx = 10;
neighbors = mysort.mea.electrodesCloserThan(channel_x, channel_y, eidx, 50);
mysort.plot.figure;
for i=1:length(channel_x)
    plot(channel_x(i), channel_y(i), '.k');
    hold on
end
plot(channel_x(eidx), channel_y(eidx), 'or');
hold on
for i=1:length(neighbors)
    eidx = neighbors(i);
    plot(channel_x(eidx), channel_y(eidx), 'ob');
end






