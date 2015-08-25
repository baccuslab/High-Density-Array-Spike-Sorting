function C = cov(f, flag);
%LDA/COV Within-groups covariance matrix of LDA object.
%   C = COV(F) returns the within-groups covariance matrix from the
%   linear discriminant analysis object F. COV(F) or COV(F, 0) returns
%   the unbiased estimate of the within-groups covariance matrix
%   regardless of the estimate obtained in F.SCALE in the function
%   LDA. If the LDA object F represents either an unbiased or
%   maximum-likelihood estimation of the within-group means and
%   covariance matrices, this will differ by a factor of (n-g)/(n-1)
%   from the value of COV(X - F.MEANS(K,:)) where X and K are the
%   original observations used in the the training of the LDA object F
%   and n and g are the number of observations and classes
%   respectively.
%
%   C = COV(F, 1) returns the maximum likelihood estimate of the
%   covariance matrix regardless of the estimate obtained in
%   F.SCALE. This should be the same as COV(X - F.MEANS(K, :), 1)
%   to within a small error.
%
%   C = COV(F, F.EST) returns either the maximum likelihood or
%   unbiased estimate depending on the original estimate of
%   F.SCALE.
%
%   If F represents a t-parameter estimation of within-groups means
%   and covariance, COV(F, 1) is determined as NU*C/(2-NU) where NU is
%   the degrees of freedom of the t-parameter estimation and C is the
%   covariance matrix of the LDA object.
%
%   C = COV(F.MEANS, ...) obtains the between-groups covariance matrix.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

error(nargchk(1, 2, nargin))
if nargin < 2
  flag = 0;
elseif isempty(flag) | (~ischar(flag) & ~isa(flag, 'double')) | ...
      length(flag) ~= 1 | ~any(flag == [0, 1, 't'])
  error('FLAG must be either 0, 1, or ''t''.')
elseif f.est == 't' & flag == 1
  flag = 't';
end

S = inv(f.scale);
C = S'*S;

if f.est == 't'
  C = f.nu*C./(f.nu-2);
end

if flag ~= f.est
  g = size(f.means, 1);
  n = sum(f.classifier.counts);
  if flag == 1
    C = C*(n - g)/n;
  else
    C = C*n/(n - g);
  end
end









