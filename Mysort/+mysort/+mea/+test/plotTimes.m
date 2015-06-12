

TIMES = abs(randn(2, 3, 5));
TIMES(1,:,:) = TIMES(1,:,:)*2.5;

figure;
ah(1) = subplot(1,2,1);
boxplot(squeeze(TIMES(1,:,:))')
ah(2) = subplot(1,2,2);
boxplot(squeeze(TIMES(2,:,:))', 'colors', 'g')
linkaxes(ah, 'xy');
set(gca, 'xtick', [1 2 3], 'xticklabel', {'a', 'b'})