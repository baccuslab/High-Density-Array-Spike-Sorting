function [c, post, V] = crossval(rule, X, k, v)
%CROSSVAL K-fold cross validation.
%   [C, POST] = CROSSVAL(RULE, X, K) returns the leave-one-out
%   cross-validation classification estimates of the observations
%   described by the feature matrix X and the class index vector K
%   given the name of the classifier RULE. RULE must be a string
%   containing the name of a CLASSIFIER object such as one of the
%   descendents of the class 'classifier'. X is an n by p matrix where
%   n is the number of observations and p is the number of
%   variables. K is a length n vector of class indeces in which each
%   element gives the group to which the corresponding row in X
%   belongs. C contains the estimated classificiations in a length n
%   class index vector and POST is a n by g matrix of posterior
%   probabilities.
%
%   [C, POST, V] = CROSSVAL(RULE, X, K, NSETS) where NSETS is a scalar
%   value no more than n perfoms NSETS-fold cross-validation with
%   randomly chosen sets of approximately equal size. The return
%   argument V is a vector of indeces indicating the cross-validation
%   group to which each row of X and K were assigned. If NSETS is a
%   vector of indeces, it is used as the cross-validition sets
%   instead. Empty sets are not permitted in NSETS.
%
%   CROSSVAL(RULE, X, K, OPTS) allows optional arguments to be passed
%   in the fields of OPTS. The only field that is used by CROSSVAL is
%   NSETS. Additionally, OPTS is passed to the CLASSIFIER constructor
%   as well as the function CLASSIFY for the appropriate CLASSIFIER
%   class so that options may be assigned to those funcitons as well.

%   Copyright (C) Michael Kiefte 1998.

%   $Log$

error(nargchk(3, 4, nargin))
h = classifier(X, k);
[n p] = size(X);
g = max(k);

if nargin < 4
  v = [];
end

if nargin > 3 & isstruct(v)
  opt = v;
  v = parseopt(opt, 'nsets');
else
  opt = [];
end

if isempty(v)
  V = 1:n;
  v = n;
elseif length(v) == 1
  if ~isa(v, 'double') | ~isreal(v) | v <= 0 | round(v) ~= v
    error('V must be a positive, non-zero integer.')
  elseif v > n
    error('V must be less than the number of observations.')
  elseif v == n
    V = 1:n;
  else
    w = v*fix(n/v);
    r = randperm(v);
    V(randperm(n)) = [reshape(repmat(1:v, w/v, 1), w, 1); ...
		      r(1:n-w)'];
  end
else
  if ~isa(v, 'double') | ~isreal(v) | prod(size(v)) ~= length(v) | ...
	any(round(v) ~= v | v <= 0 | isinf(v))
    error(['Cross-validation set indeces V must be a vector of' ...
	   ' positive, finite, non-zero integers.'])
  elseif length(v) ~= n
    error(['Length of set indeces V must be same as number of' ...
	   ' observations in X and K.'])
  elseif ~all(sum(sparse(1:n, v, 1)))
    error('Set indeces V may not have empty groups.')
  end
  
  V = v;
  v = max(V);
end
  
c = zeros(n, 1);
if nargout > 1
  post = zeros(n, g);
end

for i = 1:v
  s = warning;
  warning off
  f = feval(rule, X(V ~= i, :), k(V ~= i), opt);
  warning(s);
  if nargout == 1
    c(V == i) = classify(f, X(V == i, :), opt);
  else
    [c(V == i) post(V == i, :)] = classify(f, X(V == i, :), opt);
  end
end







