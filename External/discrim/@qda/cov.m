function c = cov(f, flag)
%QDA/COV Within-groups covariance matrices.
%   C = COV(F) returns the within-groups covariance matrices from the
%   quadratic discriminant analysis object F in the p by p by g array
%   C, where p is the number of features and g is the number of
%   classes in the original training data for F. COV(F) or COV(F, 0)
%   returns the unbiased estimates of the within-groups covariance
%   matrix regardless of the estimates obtained in F.SCALE in the
%   function QDA. If the QDA object F represents an unbiased or
%   maximum-likelihood estimation of the within-group means and
%   covariance matrices, then C(:,:,i) should be COV(X(K==i,:)) where
%   X and K represent the original training data.
%
%   C = COV(F, 1) returns the maximum likelihood estimates of the
%   covariance matrices regardless of the estimate obtained in
%   F.SCALE. Therefore, C(:,:,i) should be the same as
%   COV(X(K==i,:),1).
%
%   C = COV(F, F.EST) returns either the maximum likelihood or
%   unbiased estimates depending on the original estimate of F.SCALE.
%
%   If F represents a t-parameter estimation of within-groups means
%   and covariances, COV(F, 1) is determined as NU*C/(2-NU) where NU
%   is the degrees of freedom of the t-parameter estimation and C are
%   the covariance matrices of the QDA object.
%
%   C = COV(F.MEANS, ...) obtains the between-groups covariance matrix.
%   
%   C = COV(LDA(F), ...) returns the pooled within-groups
%   covariance matrix.

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

[g p] = size(f.means);
c = zeros(p, p, g);
nj = f.classifier.counts;
n = sum(nj);
for i = 1:g
  S = inv(f.scale(:, :, i));
  c(:,:,i) = S'*S;
end

if f.est == 't'
  c = f.nu*c./(f.nu-2);
end

if flag ~= f.est
  if flag == 1
    c = c.*repmat(shiftdim((nj-1)./nj, -1), [p p 1]);
  else
    c = c.*repmat(shiftdim(nj./(nj-1), -1), [p p 1]);
  end
end

