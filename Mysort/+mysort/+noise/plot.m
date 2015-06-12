
function C = plot(h, lims, x, idx, tau, i1, i2, col, legendstr)
set(h, 'fontsize', 14);
if ~exist('col', 'var')
    col = 'k.';
end

plot(1,1, '.w'); hold on
plot(x(i1,idx(1:end-tau)), x(i2,idx(1+tau:end)), col, 'markersize', 1);
set(h, 'xlim', lims, 'ylim', lims);
if exist('legendstr', 'var')
    legend(legendstr, 'Location', 'SouthEast');
    legend('boxoff')
end
xlabel(sprintf('x_{%d,t}', i1));
ylabel(sprintf('x_{%d,t+%d}', i2, tau));
axis equal
axis square
C = cov(x(i1,idx(1:end-tau)), x(i2,idx(1+tau:end)));
%title(sprintf('%.1f | %.1f', C(1,1), C(1,2)));
title(sprintf('cov = %.1f', C(1,2)));
set(h, 'xlim', lims, 'ylim', lims);