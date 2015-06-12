function m = confmat(c, d)
%CONFMAT Confusion matrix.
%   M = CONFMAT(C, D) returns the confusion matrix based on observed
%   class index vector C and estimated class index vector D. D may
%   contain empty classes but C may not. The result is a g by g
%   matrix of probabilities where rows correspond to actual classes
%   (in C) and columns correspond to estimated classed (in
%   D). Either C or D may actually be matrices of posterior
%   probabilities in which all rows sum to 1.
%
%   See the help for LDA for an example of how to use CONFMAT.
%
%   See also CLASSIFY, CROSSVAL, MCNEMAR.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

error(nargchk(2, 2, nargin))

if isempty(c) | ~isa(c, 'double') | ~isreal(c)
  error('C must be a real numeric double array.')
elseif prod(size(c)) ~= length(c)
  if ndims(c) ~= 2 | any(any(c < 0 | isnan(c)))
    error(['Posterior probabilities C must be a 2-d array of positive' ...
	   ' values.'])
  end
  
  [n g] = size(c);
  
  if any(abs(sum(c, 2) - 1) > g*eps)
    error('Rows of posterior probability matrix C must sum to 1.')
  end
  
  nj = full(sum(c));
  C = c;
else
  if any(round(c) ~= c | c <= 0 | isinf(c))
    error(['Class indeces C must be a vector of positive, finite,' ...
	   ' non-zero integers.'])
  end

  n = length(c);
  C = sparse(1:n, c, 1);
  nj = sum(C);
  g = length(nj);
  
  if ~all(nj)
    error('Must not have empty groups in class index C.')
  elseif length(nj) == 1
    error('Class index C must have more than 1 group.')
  end
end

if isempty(d) | ~isa(d, 'double') | ~isreal(d)
  error('D must be a real numeric double array.')
elseif prod(size(d)) ~= length(d)
  if ndims(d) ~= 2 | any(any(d < 0 | isnan(d)))
    error(['Posterior probabilities D must be a 2-d array of positive' ...
	   ' values.'])
  end
  
  if size(d, 1) ~= n
    error('Number of observations in C and D must match.')
  elseif size(d, 2) > g
    error('May not have more categories in D than C.')
  elseif size(d, 2) < g
    d(1, g) = 0;
  end
  
  if any(abs(sum(d, 2) - 1) > g*eps)
    error('Rows of posterior probability matrix D must sum to 1.')
  end
  
  D = d;
else
  if any(round(d) ~= d | d <= 0 | isinf(d))
    error(['Class indeces D must be a vector of positive, finite,' ...
	   ' non-zero integers.'])
  end

  if length(d) ~= n
    error('Number of observations in C and D must match.')
  elseif max(d) > g
    error('May not have more groups in D than C.')
  end
  
  D = sparse(1:n, d, 1);
end

m = diag(1./full(nj)) * full(C'*D); 






