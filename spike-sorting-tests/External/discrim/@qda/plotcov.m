function plotcov(f, conf)
%QDA/PLOTCOV Plot covariance ellipses for QDA object.
%   PLOTCOV(F) plots ellipsoids whose shape is determined by the
%   transformation matrices in F.SCALE and with centres determined by
%   F.MEANS.
%
%   PLOTCOV(F, CONF) where CONF is a value between 0 and 1, specifies
%   the confidence bounds for the covariance ellipsoids. CONF scales
%   the ellipsoids so that their coverage of the multivariate normal
%   distribution is approximately CONF. By default the coverage is
%   only about .4 for 2 dimensions or .2 for 3 which corresponds to no
%   resizing of the ellipsoids from their original scale matrix.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

error(nargchk(1, 2, nargin))
M = f.means;
[g p] = size(M);
if p ~= 2 & p ~= 3
  error('Can only plot QDA objects with 2 or 3 covariates.')
end

if nargin < 2
  conf = [];
end

if isempty(conf)
  crit = 1;
elseif ~isa(conf, 'double') | ~isreal(conf) | length(conf) ~= 1 | ...
      conf <= 0 | conf >= 1 | isnan(conf)
  error('CONF must be a positive non-zero scalar less than 1.')
else
  h = inline('gammainc(x/2, p/2)-q', 'x', 'p', 'q');
  gold = 2 - (3 - sqrt(5))/2;
  lim = [0 1];
  while prod(h(lim, p, conf)) > 0
    lim = [lim(2); lim(2)*gold];
  end
  crit = fzero(h, lim, optimset('display', 'off'), p, conf);
end

S = f.scale;

if p == 2
  t = linspace(0, 2*pi, 128)';
  e = [cos(t) sin(t)];
else
  col = vga;
  col = col([3:8 10:16], :);
  [x y z] = sphere(128);
  e = [x(:) y(:) z(:)];
end

np = get(gca, 'NextPlot');
for i = 1:g
  Me = e/S(:, :, i)*crit + repmat(M(i, :), size(e, 1), 1);
  if p == 2
    plot(Me(:,1), Me(:,2));
    text(M(i,1), M(i,2), '+', 'HorizontalAlignment', 'center', ...
	 'Color', 'b', 'FontSize', 24)
  else
    h = surfl(reshape(Me(:,1), 129, 129), reshape(Me(:,2), 129, 129), ...
	      reshape(Me(:,3), 129, 129));
    set(h, 'FaceColor', col(i, :), 'LineStyle', 'none');
  end
  if i == 1
    set(gca, 'NextPlot', 'add');
  end
end

set(gca, 'NextPlot', np);
grid on
xlabel('First variate')
ylabel('Second variate')
title('Quadratic discriminant analysis')
if p == 3
  zlabel('Third variate')
  view(3)
  axis vis3d
  camlight
  lighting phong
end




