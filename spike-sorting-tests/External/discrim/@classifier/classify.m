function [c, post] = classify(f, X);
%CLASSIFIER/CLASSIFY Categorise new data with CLASSIFIER object.
%   [C, POST] = CLASSIFY(F, X) classifies the rows of the n by p
%   feature matrix X given the CLASSIFIER object F, where n is the
%   number of observations or rows in X and p is the number of
%   features or variates. The estimated classes are returned in the
%   length n index vector C, while the posterior probabilities for
%   each class are given in the n by g matrix POST, where g is the
%   number of groups classifiable by F. Each row corresponds to a
%   row in X.
%
%   The generic CLASSIFIER object simply returns a column C containing
%   n copies of the class index with the highest prior probability
%   while each row of POST is a copy of the prior probability vector
%   F.prior.

%   Copyright (C) 1999 Michael Kiefte.

%   $Log$

error(nargchk(2, 2, nargin))

if isempty(X) | ~isa(X, 'double') | ~isreal(X) | ndims(X) ~= 2 | ...
      any(any(isnan(X) | isinf(X)))
  error(['Feature matrix X must be a 2-d array of real, finite' ...
	 ' values.'])
end

[n p] = size(X);
g = length(f.counts);

if p ~= size(f.range, 2)
  error(sprintf('Feature matrix X must have %d columns.', ...
		size(f.range, 2)))  
end

if nargout >= 1
  if isempty(f.prior)
    [y i] = max(f.counts);
    c = i(ones(n, 1), 1);
  elseif length(f.prior) == 1
    c = ones(n, 1);
  else
    [y i] = max(f.prior);
    c = i(ones(n, 1), 1);
  end
  
  if nargout >= 2
    if isempty(f.prior)
      post = repmat(f.counts/sum(f.counts), n, 1);
    elseif length(f.prior) == 1
      post = repmat(1/g, n, g);
    else
      post = repmat(f.prior, n, 1);
    end
  end
end
