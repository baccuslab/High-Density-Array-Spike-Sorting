

fh = figure();
ax = axes();
X = randn(100,10);
plot(X);
hold on
plot(mean(X, 2), 'k', 'linewidth', 3)
plot(mean(X, 2)+3*std(X, [], 2), ':k', 'linewidth', 3)
plot(mean(X, 2)-3*std(X, [], 2), ':k', 'linewidth', 3)

mysort.plot.legend(ax, {{'-k', 'linewidth', 3}, {':k', 'linewidth', 3}}, {'Mean', '+- 3*Std'});
mysort.plot.legend({{'-k', 'linewidth', 3}, {':k', 'linewidth', 3}}, {'Mean', '+- 3*Std'});
mysort.plot.legend(ax, {{'-k', 'linewidth', 3}, {':k', 'linewidth', 3}}, {'Mean', '+- 3*Std'}, 'NorthEastOutside');
mysort.plot.legend({{'-k', 'linewidth', 3}, {':k', 'linewidth', 3}}, {'Mean', '+- 3*Std'}, 'NorthEastOutside');



