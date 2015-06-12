function plotdr(f, varargin)
%PLOTDR Plot decision regions for classifier object.
%   PLOTDR(F, ...) plots the decision boundaries of maximum posterior
%   likelihood for different classes where F is an object descended
%   from CLASSIFIER. Additional arguments are passed to the
%   approprate CLASSIFY function.
%
%   PLOTDR can plot decision regions for classifiers with two variates
%   only.
%
%   See the help for LOGDA for an example of how to use PLOTDR.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

error(nargchk(1, inf, nargin))

p = f.nvar;
if p ~= 2
  error('Can only plot decision region for 2 covariates.')  
end
g = f.nclass;

np = get(gca, 'NextPlot');
if strcmp(np, 'add')
  a = axis;
else
  a = f.range;
  h = plot(a(:,1), a(:,2));
  a = axis;
  delete(h);
end
nint = 128;

x = linspace(a(1), a(2), nint);
y = linspace(a(3), a(4), nint);
[X Y] = meshgrid(x, y);

[c post] = classify(f, [X(:) Y(:)], varargin{:});

for k = 1:g
  s = reshape(max(post(:, [1:k-1 k+1:g]), [], 2) - post(:,k), ...
	      repmat(nint, 1, p));
  contour(x, y, s, [0 0], 'b');
  if k == 1
    set(gca, 'NextPlot', 'add');
  end
end
set(gca, 'NextPlot', np);

grid on
xlabel('First variate')
ylabel('Second variate')
title('Decision regions')
axis(a)
