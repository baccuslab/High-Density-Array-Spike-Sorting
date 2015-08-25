function plotobs(X, k)
%PLOTOBS Scatter plot of observation classes
%   PLOTOBS(X, K) where X is a 2 or 3 column matrix of features and K
%   is a vector of class indeces produces a scatter plot of
%   observations by class. The default labels are 'a' through 'z'. X
%   must have at least 2 and no more than 3 columns or features.
%
%   See help for LOGDA for an example of how to use PLOTOBS.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

error(nargchk(2, 2, nargin))

if isempty(X) | ~isa(X, 'double') | ~isreal(X) | ndims(X) ~= 2 | ...
      any(any(isnan(X) | isinf(X)))
  error(['Feature matrix X must be a 2-d array of real, finite' ...
	 ' values.'])
end

[n p] = size(X);

if p == 1 | p > 3
  error('Can only plot observations with 2 or 3 variables.')
end

if isempty(k) | ~isa(k, 'double') | ~isreal(k) | ...
      prod(size(k)) ~= length(k) | ...
      any(round(k) ~= k | k <= 0 | isinf(k))
  error(sprintf(['Observed class indeces K must be a vector of' ...
		 ' positive, finite,\nnon-zero integers.']))  
elseif length(k) ~= n
  error(['Class index K must have same number of observations' ...
	 ' as\nobserved class matrix K.'])  
end

g = max(k);
if g > 26
  error('Too many classes.')
end

switch p
 case 2
  h = plot(X(:,1), X(:,2));
  a = axis;
  delete(h);
  for i = 1:g
    text(X(k == i, 1), X(k == i, 2), char('a' + i - 1), ...
	 'HorizontalAlignment', 'center')
  end
  axis(a);
 case 3
  h = plot3(X(:,1), X(:,2), X(:,3));
  a = axis;
  delete(h);
  for i = 1:g
    text(X(k == i, 1), X(k == i, 2), X(k == i, 3), ...
	 char('a' + i - 1), 'HorizontalAlignment', 'center')
  end
  axis(a)
  axis vis3d
  zlabel('Third variate')
end

xlabel('First variate')
ylabel('Second variate')
grid on







