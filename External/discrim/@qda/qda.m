function [f, c, post] = qda(X, k, prior, est, nu)
%QDA Quadratic Descriminant Analysis.
%   F = QDA(X, K, PRIOR) returns a quadratic discriminant analysis
%   object F based on the feature matrix X, class indeces in K and the
%   prior probabilities in PRIOR where PRIOR is optional. See the help
%   for QDA's parent object CLASSIFIER for information on the input
%   arguments X, K and PRIOR.
%
%   In addition to the fields defined by the CLASSIFIER class, F
%   contains the following fields:
%
%   MEANS: a g by p matrix where g is the number of classes and p is
%   the number of features or variates. Each row gives the mean vector
%   for each class. 
%
%   SCALE: the p by p by g numeric array in which each p by p matrix
%   is the scale matrix that transforms the observed within-groups
%   covariance for the corresponding class to identity. Therefore
%   F.SCALE(:,:,i)=INV(CHOL(COVX(K==i,:)),1) for maximum-likelihood
%   estimates (see below) or INV(CHOL(COV(X(K==1,:)))) for ubiased
%   estimates.
%
%   LDET: the length g vector which gives the log determinants for
%   each covariance matrix.
%
%   EST: either 0, 1, or 't' representing unbiased, maximum likelihood
%   or t-parameter estimation respectively as explained below.
%
%   NU: This field is only present if EST is 't'. NU gives the degrees
%   of freedom for the t-parameter estimation as explained in the next
%   paragraph.
%
%   QDA(X, K, PRIOR, EST, NU) where EST is one of 'unbiased', 'ml', or
%   't', uses either bias-corrected, maximum likelihood or t-parameter
%   estimation respectively. For t-parameter estimation, an additional
%   argument, NU, gives the degrees of freedom for the estimator (the
%   default is 5 if not given). The default estimator is unbiased
%   estimation (which corresponds to the default for the functions STD
%   and COV). Unbiased estimation bias corrects the estimate for the
%   within-groups covariance matrix by a factor of 1/(n(i)-1) where
%   n(i) is the number of observations in class i (as returned by
%   F.COUNTS). For maximum likelihood estimation, no correction is
%   made. For t-parameter estimation, the means and scale matrix are
%   estimated by an iterative weighted algorithm. When specifying EST,
%   only the first few disambiguating letters need be given: i.e.,
%   'u', 'm' or 't'.
%
%   QDA(X, K, EST) is equivalent to QDA(X, K, [], EST).
%
%   QDA(X, K, OPTS) allows optional arguments to be passed in the
%   fields of the structure OPTS. Fields that are used by QDA are
%   PRIOR, EST and NU.
%
%   [F, C, POST] = QDA(X, K, ...) additionally performs leave-one-out
%   cross-validation on the data in X. C is a length n index vector of
%   estimated class memberships similar to K corresponding to the
%   matrix of features X. POST is an n by g matrix of posterior
%   probabilities. Leave-one-out cross-validation is only defined for
%   methods 'ml' and 'unbiased'. C and POST will not necessarily
%   correspond to the output of CROSSVAL(X, K, 'qda', ...) because in
%   the latter, the prior probabilities are not fixed between
%   cross-validation estimates unless this is done so explicitly in
%   the option struct passed to CROSSVAL.
%
%   See also CLASSIFIER, LDA, LOGDA, SOFTMAX, COV, CROSSVAL.
%
%   Example:
%   %generate artificial data with 4 classes and 3 variates
%   r = randn(3, 3, 4);
%   for i = 1:4
%      % generate random covariance matrices for each class
%      C(:,:,i) = r(:,:,i)'*r(:,:,i);
%   end
%   M = randn(4, 3)*2;       % random means
%   k = ceil(rand(400, 1)*4);  % random classes
%   X = randn(400, 3);
%   for i = 1:4
%      X(k==i,:) = X(k==i,:)*chol(C(:,:,i)) + M(k(k==i), :);
%   end
%   f = qda(X, k); disp(f)
%   cov(f), plotcov(f)
%   plotcov(shrink(f, 1))
%   g = lda(X, k);
%   [m alpha] = mcnemar(k, f(X), g(X))
%
%   References:
%   B. D. Ripley (1996) Pattern Classification and Neural
%   Networks. Cambridge.

%   Copyright (c) 1999 Michael Kiefte.  
%   Additionally based on algorithm presented in S-Plus code written
%   by Ripley and Venables.

%   $Log$

error(nargchk(2, 5, nargin))

if nargin > 2 & isstruct(prior)
  if nargin > 3
    error(sprintf(['Cannot have arguments following option struct:\n' ...
		   '%s'], nargin(3, 3, 4)))
  end
  [prior est nu] = parseopt(prior, 'prior', 'est', 'nu');
elseif nargin < 5
  nu = [];
  if nargin < 4
    est = [];
    if nargin < 3
      prior = [];
    end
  end
end

if ischar(prior)
  nu = est;
  est = prior;
  prior = [];
end

[h G] = classifier(X, k, prior);
[n p] = size(X);
nj = h.counts;
g = length(nj);
prior = h.prior;

if nargout > 1
  cv = 1;
else
  cv = 0;
end

if isempty(est)
  est = 0;
elseif ~ischar(est) | length(est) ~= size(est, 2) | ...
      size(est, 1) ~= 1
  error('EST must be a string.')
else
  t = find(strncmp(est, {'unbiased', 'ml', 't'}, length(est)));
  if isempty(t)
    error('EST must be one of ''unbiased'', ''ml'', or ''t''.')
  end
  switch t
   case 1
    est = 0;
   case 2
    est = 1;
   otherwise
    est = 't';
  end
end

if est == 't'
  if isempty(nu)
    nu = 5;
  elseif ~isa(nu, 'double') | length(nu) ~= 1 | round(nu) ~= nu | ...
	nu < 3 | isinf(nu)
    error(['Degrees of freedom NU must be a finite, integer scalar' ...
	   ' greater than 2.'])
  elseif cv
    error('Cannot perform cross-validation with t-estimator.')
  end
elseif ~isempty(nu)  
  error('May specify degrees of freedom NU only with t-estimator.')  
end

M = sparse(1:g, 1:g, 1./nj')*G'*X;
S = zeros(p, p, g);
ldet = zeros(1, g);

for i = 1:g
  switch est
   case {0, 1}
    r = qr((X(k == i,:) - repmat(M(i,:), nj(i), 1)) ...
	       /sqrt(nj(i) - (1-est)));
   otherwise
    w = ones(nj(i), 1);
    Xk = X(k == i,:);
    c = (nu+p)/(nj(i)*nu);
    while 1
      wold = w;
      Xc = Xk - repmat(M(i,:), nj(i), 1);
      r = triu(qr(repmat(sqrt(w*c), 1, p).*Xc));
      w = 1./(1+(Xc/r(1:p,:)).^2*repmat(1/nu, p, 1));
      M(i,:) = w'*Xk/sum(w);
      if max(abs(w-wold)) < max(w)*nj(i)*eps 
	break
      end
    end
  end
  S(:,:,i) = inv(triu(r(1:p,:)));
  ldet(i) = 2*sum(log(abs(diag(r))));
end

if cv
  lc = ldet(ones(n, 1), :);
  D = zeros(n, g);
  for i = 1:g
    D(:,i) = sum(((X - M(i(ones(n, 1)), :)) * S(:,:,i)).^2, 2);
  end
  K = 1-est;
  nc = nj(k)';
  idx = (k-1)*n+(1:n)';
  lc(idx) = lc(idx) + p*log((nc - K)./(nc - 1 - K)) + ...
	    log(1 - nc./((nc - 1).*(nc - K)).*D(idx));
  D(idx) = D(idx) .* (nc.^2.*(nc - 1 - K)) ./ ...
	   ((nc - 1).^2.*(nc - K)) ./ ...
	   (1 - nc./((nc - 1).*(nc - K)).*D(idx));
  D = (D + lc)/2 - repmat(log(prior), n, 1);
  [y c] = min(D, [], 2);
  if nargout > 2
    D = exp(y(:, ones(1, g)) - D);
    post = D./repmat(sum(D, 2), 1, g);
  end
end

f = class(struct('means', M, 'scale', S, 'ldet', ldet, 'est', est, ...
		 'nu', nu), 'qda', h);

