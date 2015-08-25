function [m, alpha, na, nb] = mcnemar(k, a, b)
%MCNEMAR McNemar's test for comparison of error rates.
%   [M, ALPHA, NA, NB] = MCNEMAR(K, A, B) returns the McNemar's test
%   statistic for the comparison of error rates in A versus B given
%   the actual observed categories in K. K, A, and B are length n
%   vectors of class indeces. K may not have empty classes. NA is the
%   number of errors in A that are not made in B and NB is the number
%   of errors in B that are not made in A. M is then 
%
%   (ABS(NA-NB) - 1)^2/(NA+NB) 
%
%   which can then be compared to a chi squared distribution on 1
%   degree of freedom as a test for the improvement in correct
%   classification in A versus B. ALPHA is the probability of
%   observing a value of NA or less given the null hypothesis binomial
%   distribution of B(NA+NB,1/2) and serves as a test for the
%   improvement of the estimation in A over B.
%
%   Example:
%   [m, alpha] = mcnemar(k, a, b);
%   1 - gammainc(m/2, .5)    %chi squared test of m
%
%   References:
%   J. L. Fleiss (1981) Statistical Methods for Rates and
%   Proportions. Second Edition. Wiley.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

error(nargchk(3, 3, nargin))

if isempty(k) | ~isa(k, 'double') | ~isreal(k) | ...
      prod(size(k)) ~= length(k) | ...
      any(round(k) ~= k | k <= 0 | isinf(k))
  error(sprintf(['Observed class indeces K must be a vector of' ...
		 ' positive, finite,\nnon-zero integers.']))  
end

n = length(k);
if ~all(sum(sparse(1:n, k, 1)))
  error('Must not have empty groups in class index K.')
end
g = max(k);
  
if isempty(a) | ~isa(a, 'double') | ~isreal(a) | ...
      prod(size(a)) ~= length(a) | ...
      any(round(a) ~= a | a <= 0 | isinf(a))
  error(sprintf(['Predicted class indeces A must be a vector of' ...
		 ' positive, finite,\nnon-zero integers.']))  
elseif length(a) ~= n
  error(['Predicted class index A must have same number of observations' ...
	 ' as\nobserved class matrix K.'])  
elseif max(a) > g
    error('May not have more groups in predicted class index A than in K.')    
end

if isempty(b) | ~isa(b, 'double') | ~isreal(b) | ...
      prod(size(b)) ~= length(b) | ...
      any(round(b) ~= b | b <= 0 | isinf(b))
  error(sprintf(['Predicted class indeces B must be a vector of' ...
		 ' positive, finite,\nnon-zero integers.']))  
elseif length(b) ~= n
  error(['Predicted class index A must have same number of observations' ...
	 ' as\nobserved class matrix K.'])  
elseif max(b) > g
    error('May not have more groups in predicted class index A than in K.')    
end

na = sum(a ~= k & b == k);
nb = sum(b ~= k & a == k);
m = (abs(na - nb) - 1)^2/(na + nb);

if nargout >= 2
  n = na+nb;
  s = 0:min(na, nb-1);
  alpha = sum(exp(gammaln(n+1) - gammaln(s+1) - ...
		  gammaln(n - s + 1))) * .5^n;
  if na >= nb
    alpha = 1 - alpha;
  end
end






