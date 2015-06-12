function [f, c, post] = lda(X, k, prior, est, nu)
%LDA Linear Discriminant Analysis.
%   F = LDA(X, K, PRIOR) returns a linear discriminant analysis object
%   F based on the feature matrix X, class indeces in K and the prior
%   probabilities in PRIOR where PRIOR is optional. See the help for
%   LDA's parent object CLASSIFIER for information on the input
%   arguments X, K and PRIOR.
%
%   In addition to the fields defined by the CLASSIFIER class, F
%   contains the following fields:
%
%   MEANS: a g by p matrix where g is the number of classes and p is
%   the number of features or variates. Each row gives the mean vector
%   for each class. 
%
%   SCALE: the p by p matrix which transforms the observed
%   within-groups covariance to identity. This is equivalent to
%   INV(CHOL(COV(X - F.MEANS(K, :), 1))) for maximum-likelihood
%   estimates (see below). For unbiased estimates,
%   INV(CHOL(COV(F.MEANS(K, :)))) will differ from F.SCALE by a factor
%   of SQRT((n-1)/(n-g)) because the two normalize on a different
%   number of estimated means.
%
%   EST: either 0, 1, or 't' representing unbiased, maximum likelihood
%   or t-parameter estimation respectively as explained below.
%
%   NU: This field is only present if EST is 't'. NU gives the degrees
%   of freedom for the t-parameter estimation as explained in the next
%   paragraph.
%
%   LDA(X, K, PRIOR, EST, NU) where EST is one of 'unbiased', 'ml', or
%   't', uses either bias-corrected, maximum likelihood or t-parameter
%   estimation respectively. For t-parameter estimation, an additional
%   argument, NU, gives the degrees of freedom for the estimator (the
%   default is 5 if not given). The default estimator is unbiased
%   estimation (which corresponds to the default for the MATLAB
%   functions STD and COV). Unbiased estimation bias corrects the
%   estimate of the within-groups covariance matrix by a factor of
%   1/(n-g). For maximum likelihood estimation, no correction is
%   made. For t-parameter estimation, the means and scale matrix are
%   estimated by an iterative weighted algorithm. When specifying EST,
%   only the first few disambiguating letters need be given: i.e.,
%   'u', 'm' or 't'.
%
%   LDA(X, K, EST) is equivalent to LDA(X, K, [], EST).
%
%   LDA(X, K, OPTS) allows optional arguments to be passed in the
%   fields of the structure OPTS. Fields that are used by LDA are
%   PRIOR, EST and NU.
%
%   [F, C, POST] = LDA(X, K, ...) additionally performs leave-one-out
%   cross-validation on the data in X. C is a length n index vector of
%   estimated class memberships similar to K corresponding to the
%   matrix of features X. POST is an n by g matrix of posterior
%   probabilities. Leave-one-out cross-validation is only defined for
%   methods 'ml' and 'unbiased'. C and POST will not necessarily
%   correspond to the output of CROSSVAL(X, K, 'lda', ...) because in
%   the latter, the prior probabilities are not fixed between
%   cross-validation estimates unless this is done so explicitly in
%   the option struct passed to CROSSVAL.
%
%   F = LDA(G) where G is an object of class QDA returns an LDA object
%   based on G.
%
%   See also CLASSIFIER, QDA, LOGDA, SOFTMAX, COV, CROSSVAL.
%
%   Example:
%   %generate artificial data with 4 classes and 3 variates
%   r = randn(3);
%   C = r'*r;     % random positive definite symmetric matrix
%   M = randn(4, 3)*2;    % random means
%   k = ceil(rand(400, 1)*4);    % random class indeces
%   X = randn(400, 3)*chol(C) + M(k, :);  
%   f = lda(X, k); disp(f)
%   [lambda ratio] = cvar(f)   % canonical variates
%   cov(f), plotcov(f)
%   plotcov(shrink(f, .5))
%   [c post] = classify(f, X);
%   confmat(k, c)
%   confmat(k, post)
%
%   References: 
%   B. D. Ripley (1996) Pattern Classification and Neural
%   Networks. Cambridge.

%   Copyright (c) 1999 Michael Kiefte.  

%   $Log$

if isa(X, 'qda')
  error(nargchk(1, 1, nargin))
  g = shrink(X, 1);
  if g.est == 't'
    warning(['Not an exact conversion for t-estimator QDA' ...
	     ' objects.'])
    nu = g.nu;
  else
    nu = [];
  end

  f = class(struct('means', g.means, 'scale', g.scale(:,:,1), ...
		   'est', g.est, 'nu', nu), 'lda', g.classifier);
  return
end

error(nargchk(2, 5, nargin))

if nargin > 2 & isstruct(prior)
  if nargin > 3
    error(sprintf(['Cannot have arguments following option struct:\n' ...
		   '%s'], nargchk(3, 3, 4)))    
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

M = sparse(1:g, 1:g, 1./nj)*G'*X;
Xc = X - M(k, :);

S = std(Xc);
if any(S < n*max(S)*eps)
  error(sprintf(['Column %d in feature matrix X is constant within' ...
		 ' groups.'], min(find(S < n*max(S)*eps)))) 
end
S = diag(S);

switch est
 case {0, 1}
  [u s v] = svd(Xc*S/sqrt(n - g*(1-est)), 0);
  r = sum(diag(s) > n*s(1)*eps);
  if (r < p)
    warning(sprintf(['Nullity of within-groups covariance matrix is' ...
		     ' %d.'], p - r))
    v = v(:,1:r);
    s = s(1:r,1:r);
  end  
  S = S*inv(triu(qr(s*v')));

  if cv
    Xs = X*S;
    Ms = M*S;
    XM = Xs - Ms(k, :);
    nc = nj(k)';
    K = g*(1-est);
    c = (n - K - 1)/(n - K);
    D = repmat(sum(Xs.^2, 2), 1, g) - 2*Xs*Ms' + repmat(sum(Ms'.^2), n, 1);
    Dc = D((k-1)*n+(1:n)');
    cc = (n - K)*(nc - 1)./nc;
    D = c * (D + (repmat(sum(Xs .* XM, 2), 1, g) - XM*Ms').^2 ./ ...
	      repmat(cc - Dc, 1, g));
    D((k-1)*n+(1:n)') = Dc * c .* (nc./(nc - 1)).^2 ./ (1 - Dc./cc);
    D = D/2 - repmat(log(prior), n, 1);
    [y c] = min(D, [], 2);
    if nargout > 2
      D = exp(y(:, ones(1, g)) - D);
      post = D./repmat(sum(D, 2), 1, g);
    end
  end
 otherwise
  w = ones(n,1);
  c = (nu+p)/(n*nu);
  sing = 0;
  while 1
    wold = w;
    [u s v] = svd(repmat(sqrt(w*c), 1, p).*Xc*S, 0);
    r = sum(diag(s) > n*s(1)*eps);
    if r < p
      if ~sing
	warning(sprintf(['Nullity of within-groups covariance matrix is' ...
			 ' %d.'], p - r))	
	sing = 1;
      end
      v = v(:, 1:r);
      s = s(1:r, 1:r);
    end
    w = 1./(1+(Xc*S*v/s).^2*repmat(1/nu, p, 1));
    M = G'*(w(:,ones(1,p)).*X)./repmat(G'*w, 1, p);
    if max(abs(w - wold)) < max(w)*n*eps 
      break
    end
    Xc = X - M(k, :);
  end
  S = S*inv(triu(qr(s*v')));
end

f = class(struct('means', M, 'scale', S, 'est', est, 'nu', nu), ...
	  'lda', h);






