function [f, G, w] = classifier(X, k, prior)
%CLASSIFIER Generic discriminant analysis object.
%   F = CLASSIFIER(X, K, PRIOR) returns a generic discriminant
%   analysis object based on the feature matrix X, class indeces in K
%   and the prior probabilities in PRIOR. X must be an n by p real
%   matrix where n is the number of rows or observations and p is the
%   number of features or variates. It is an error to give an X in
%   which any of the columns are constant. K must be a vector of class
%   indeces of length n where each element gives the class to which
%   the corresponding observation in X belongs. It is an error to
%   supply a K that has empty classes: i.e. the number of unique
%   elements in K must also be the maximum value. It is also an error
%   to supply a K with only one class. The argument PRIOR is optional.
%
%   If given, PRIOR must be a length g vector of prior probabilities
%   where g is the number of classes in K, or the scalar value 1. If
%   PRIOR is a vector, it must sum to 1. Otherwise, if PRIOR is the
%   scalar value 1, equal priors are assigned to each of the g
%   classes. If not given or empty, the default is to estimate prior
%   probabilities based on the counts of each class in K.
%
%   The returned object, as well as those returned by objects decended
%   from the CLASSIFIER object class, may be used as a function
%   similar to CLASSIFY. However, a generic CLASSIFIER object does
%   nothing interesting on its own. It is the parent class to all
%   other classifier objects.
%
%   [F, G] = CLASSIFIER(X, K, PRIOR) also returns the incidence matrix
%   G based on the class indeces in K. G is an n by g sparse matrix in
%   which each row contains only one non-zero value in the column
%   indicating to which class the corresponding observation in X
%   belongs. This is usefull for some descendent classifier
%   constructors that require this matrix such as LDA.
%
%   [F, G, W] = CLASSIFER(X, G) where G is an n by g matrix of either
%   posterior probabilities or multinomial counts, attempts to model
%   these instead of the class vector K. If G is a matrix of posterior
%   probabilities, each row must sum to 1, otherwise, G must contain
%   positive integers indicating the number of class counts for each
%   group and for each corresponding predictor in X. This variation is
%   only applicable to some CLASSIFIER descendents such as LOGDA and
%   it is an error to supply the third argument PRIOR. W is the sum of
%   the rows of G and is usefull for some descendent classifier
%   constructors. The returned G is normalised such that the rows sum
%   to 1.
%
%   The fields of F are fully accessable, but may not be assigned
%   to. F contains the following fields:
%
%   PRIOR: a length g vector giving the prior probabilities for each
%   class in K. If a matrix of posterior probabilities or counts are
%   given, PRIOR is estimated from the normalised sum over rows.
%
%   COUNTS: a length g vector of total counts for each class in K. If
%   a matrix of posterior probabilities is given, this information is
%   unavailable. In this case COUNTS is a scalar integer giving the
%   total number of observations in X and G.
%
%   RANGE: a 2 by p matrix which gives the minimum and maximum for
%   each variate in the observation matrix X. The first row contains
%   the minima and the second contains the maxima.
%
%   NVAR: the number of variates or features in the observation
%   matrix X.
%
%   NCLASS: the number of groups in the class vector K or matrix G.
%
%   See also: LDA, QDA, LOGDA, SOFTMAX.

%   Copyright (C) 1999 Michael Kiefte.

%   $Log$

error(nargchk(2, 3, nargin))

if isempty(X) | ~isa(X, 'double') | ~isreal(X) | ndims(X) ~= 2 | ...
      any(any(isnan(X) | isinf(X)))
  error(['Feature matrix X must be a 2-d array of real, finite' ...
	 ' values.'])
end

[n p] = size(X);
range = [min(X); max(X)];
if ~all(diff(range))
  error(sprintf('Column %d in feature matrix is constant.', ...
		min(find(~all(diff(range))))))  
end

if isempty(k) | ~isa(k, 'double') | ~isreal(k)
  error('Class index K must contain all real, numeric values.')
elseif prod(size(k)) ~= length(k)
  if ndims(k) ~= 2 | any(any(k < 0 | isnan(k)))
    error(['Incidence matrix G must be a 2-d array of positive' ...
	   ' values.'])
  elseif size(k, 1) ~= n
    error(['G must have same number of observations as feature' ...
	   ' matrix X.'])    
  elseif nargin > 2
    error(sprintf(['May not give prior probabilities when second' ...
		   ' argument is\na matrix of posterior probabilities' ...
		   ' or an incidence matrix.']))
  end
  
  g = size(k, 2);
  
  if all(all(round(k) == k))
    nj = sum(k);
    if ~all(nj)
      error(sprintf('Column %d of G has no observations.', ...
		    min(find(~nj))))      
    elseif any(isinf(nj))
      error('Infinities in incidence matrix G')
    elseif ~all(any(k, 2))
      error(sprintf('Row %d of g has no responses.', ...
		    min(find(~all(any(k, 2))))))
    end
    nj = full(nj);
    prior = nj/sum(nj);
  elseif any(abs(sum(k, 2) - 1) > g*eps)
    error('Rows of posterior probability matrix G must sum to 1.')
  else
    prior = full(sum(k))/n;
    nj = n;
  end
  w = full(sum(k, 2));
  G = sparse(1:n, 1:n, 1./w) * k;
elseif any(round(k) ~= k | k <= 0 | isinf(k))
  error(['Class indeces K must contain all finite, positive, non-zero' ...
	 ' integers.'])
elseif length(k) ~= n
  error(['Class index K must have same number of observations as' ...
	 ' feature matrix X.'])
else
  G = sparse(1:n, k, 1);
  nj = sum(G);
  g = length(nj);
  
  if ~all(nj)
    error('Must not have empty groups in class index K.')
  elseif length(nj) == 1
    error('Class index K must have more than 1 group.')
  else
    nj = full(nj);
  end
  
  if nargin < 3
    prior = [];
  end
  
  if ~isempty(prior)
    if ~isa(prior, 'double') | ~isreal(prior) | ...
	  prod(size(prior)) ~= length(prior) | ...
	  any(prior <= 0 | isnan(prior))
      error(['Prior probabilities PRIOR must be a vector of positive' ...
	     ' values.'])
    elseif abs(sum(prior) - 1) > length(prior)*eps
      error(['Prior probabilities PRIOR must sum to 1 (or be the scalar' ...
	     ' value 1).'])
    elseif length(prior) > 1 & length(prior) ~= g
      error(['Prior probabilities PRIOR must have same number of' ...
	     ' groups as class index vector K.'])    
    end
  end
end

f = class(struct('prior', prior(:)', 'counts', nj, 'range', range), ...
	  'classifier');




