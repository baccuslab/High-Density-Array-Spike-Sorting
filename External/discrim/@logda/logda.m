function [f, iter, dev, hess] = logda(X, k, prior, maxit, est)
%LOGDA Logistic Discriminant Analysis.
%   F = LOGDA(X, K, PRIOR) returns a logistic discriminant analysis
%   object F based on the feature matrix X, class indeces in K and the
%   prior probabilities in PRIOR where PRIOR is optional. See the help
%   for LOGDA's parent object CLASSIFIER for information on the input
%   arguments X, K and PRIOR.
%
%   In addition to the fields defined by the CLASSIFIER class, F
%   contains the following field:
%
%   COEFS: a g-1 by p+1 matrix where g is the number of groups in the
%   class index K and p is the number of features or columns in
%   X. These are the coefficients from the multiple logistic
%   regression of the feature matrix on the n by g indicator matrix G
%   of K. The parameters for the first class (indexed by 1) in K is
%   assumed to have all 0 coefficients and is not represented in the
%   array. The ceofficients for the remaining classes are given in the
%   rows. The first column represents the intercept or bias term and
%   the remaining columns represent the p variables or features.
%
%   LOGDA(X, K, PRIOR, MAXITER) where MAXITER is a positive integer
%   aborts the algorithm after that many iterations. The default value
%   depends on the algorithm used (described below).
%
%   LOGDA(X, K, PRIOR, MAXITER, EST) where EST is one of 'glm' or
%   'slp' uses either a generalised linear model or a single layer
%   perceptron to minimise the residual deviance of the model. For a
%   generalised linear model, the residual deviance is minimised using
%   an iterative weighted non-linear regression using either a
%   binomial link and variance function for a two class problem or the
%   conditional poisson link and variance function for more than
%   two. The single layer perceptron uses a variable metric conjugate
%   gradient descent algorithm in which the conditional posterior
%   probabilities are used to determine the partial derivatives (see
%   SOFTMAX). For two categories, the default EST is 'glm' while for
%   more than two classes 'slp' is used unless the Hessian matrix is
%   specifically requested, in which case the default EST is
%   'slp'. Usually 'slp' is faster for large problems with more than
%   two problems, but 'glm' may be optimal for small data sets. For
%   data sets with only two classes, 'glm' is always faster. When
%   specifying EST, only the first few disambiguating letters need be
%   given: i.e., 'g' or 's'.
%
%   Either MAXITER or EST may be omitted in which case default
%   values are assigned.
%
%   LOGDA(X, K, OPTS) allows optional arguments to be passed in the
%   fields of the structure OPTS. Fields that are used by LDA are
%   PRIOR, MAXITER and EST.
%
%   [F, ITER, DEV, HESS] = LOGDA(X, K, ...) additionally returns the
%   number of iterations required by the algorithm before convergence
%   in ITER, the residual deviance for the fit in DEV, and the Hessian
%   matrix of the coefficients in HESS. HESS is a square matrix with
%   (g-1)*(p+1) rows and columns. The coefficients are ordered with
%   group categories varying fastest and with the first variate
%   representing the bias or offset term (i.e., as if the matrix
%   F.COEFS were vectorised F.COEFS(:)). The eigenvalues of HESS
%   should be all positive indicating convergance to a global
%   minimum. In order to return HESS, EST must be 'slp' if given.
%
%   LOGDA(X, G, MAXIT, EST) where G is a p by g matrix of posterior
%   probabilities or counts, models this instead of absolute class
%   memberships. If G represents counts, all of its values must be
%   positive integers. Otherwise the rows of G represent posterior
%   probabilities and must all sum to 1. It is an error to give the
%   argument PRIOR in this case. If G represents posterior
%   probabilities, F.PRIOR will be calculated as the normalised sum of
%   the columns of G and F.COUNTS will be a scalar value representing
%   the number of observations. Otherwise, F.COUNTS will be the sum of
%   the columns and F.PRIOR will represent the observed prior
%   distribution.
%
%   See also CLASSIFER, LDA, QDA, SOFTMAX.
%
%   Example:
%   %generate artificial data with 3 classes and 2 variates
%   r = randn(2);
%   C = r'*r;     % random positive definite symmetric matrix
%   M = randn(3, 2)*2;    % random means
%   k = ceil(rand(400, 1)*3);    % random class indeces
%   X = randn(400, 2)*chol(C) + M(k, :);  
%   f = logda(X, k); disp(f)
%   g = lda(X, k);
%   plotobs(X, k); hold on, plotdr(f)
%   plotdr(g), plotcov(g)
%
%   References:
%   P McCullagh and J. A. Nelder (1989) Generalized Linear
%   Models. Second Edition. Chapman & Hall.
%   B. D. Ripley (1996) Pattern Classification and Neural
%   Networks. Cambridge.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

error(nargchk(2, 5, nargin))

if nargin > 2 & isstruct(prior)
  if nargin > 3
    error(sprintf(['Cannot have arguments following option struct:\n' ...
		   '%s'], nargchk(3, 3, 4)))    
  end
  [prior maxit est] = parseopt(prior, 'prior', 'maxiter', 'est');
elseif nargin < 5
  est = [];
  if nargin < 4
    maxit = [];
    if nargin < 3
      prior = [];
    end
  end
end

if nargin >= 3 & isa(prior, 'double') & length(prior) == 1 & ...
	   prior ~= 1
  est = maxit;
  maxit = prior;
  prior = [];
elseif nargin == 3 & ischar(prior)
  est = prior;
  prior = [];
elseif nargin == 4 & ischar(maxit)
  est = maxit;
  maxit = [];
end

[n p] = size(X);

if prod(size(k)) ~= length(k)
  if ~isempty(prior)
    error(sprintf(['Cannot give prior probabilities PRIOR with' ...
		   ' incidence matrix\nor posterior probabilities G:\n' ...
		   '%s'], nargchk(2, 3, 4)))
  end
  
  [h G w] = classifier(X, k);
  g = size(G, 2);
  logG = G;
  logG(find(G)) = log(G(find(G)));
else
  [h G] = classifier(X, k, prior);
  nj = h.counts;
  g = length(nj);
  w = (nj./(n*h.prior))';
  w = w(k);
  logG = 0;
end

if isempty(est)
  if nargout == 4
    est = 1;
  elseif g == 2
    est = 0;
  else
    est = 1;
  end
elseif ~ischar(est) | ndims(est) ~= 2 | length(est) ~= size(est, 2) ...
      | size(est, 1) ~= 1
  error('EST must be a string.')
else
  est = find(strncmp(est, {'glm', 'slp'}, length(est)));
  if isempty(est)
    error('EST must be one of ''glm'' or ''slp''.')
  end
  est = est - 1;
  
  if nargout == 4 & ~est
    error('Must use single layer perceptron to output Hessian.')
  end
end

if isempty(maxit)
  if est
    maxit = 100;
  else 
    maxit = 10;
  end
elseif ~isa(maxit, 'double') | length(maxit) ~= 1 | ~isreal(maxit) | ...
      round(maxit) ~= maxit | maxit < 0
  error('Maximum number of iterations MAXITER must be a positive integer.')
end

range = h.range;
if est
  X = (X - range(ones(n,1),:)) * diag(1./diff(range));
end

trace = ~strcmp(warning, 'off');

if est
  col = sparse(n+1:g*n, repmat(1:g-1, n, 1), 1);
  U = [col [col(:, repmat((1:g-1)', 1, p)) .* ...
	    repmat(X(:,repmat(1:p, g-1, 1)), g, 1)]];
  delta = w(:, ones(1, g-1)) .* (1/g - G(:,2:end));
  grad = sum([delta, delta(:,repmat((1:g-1)', 1, p)) .* ...
	      X(:, repmat(1:p, g-1, 1))])';
  H = eye((g-1)*(p+1));
  Dold = sum(w' * (G .* (logG + log(g)) - G + 1/g));
  betaold = zeros((g-1)*(p+1), 1);
elseif g == 2
  U = [ones(n, 1) X];
  mu = (G(:,2) + .5)/2;
  eta = log(mu./(1 - mu));
  eeta = exp(eta);
  Dold = sum(w'*(G.*(logG + log(2))));
else
  col = sparse(n+1:g*n, repmat(1:g-1, n, 1), 1);
  U = [col [col(:, repmat((1:g-1)', 1, p)) .* ...
	    repmat(X(:,repmat(1:p, g-1, 1)), g, 1)], ...
      sparse(1:g*n, repmat((1:n)', 1, g), 1)];
  mu = full(G + .1);
  eta = log(mu);
  Dold = sum(w'*(G.*(logG + log(g)) - G + 1/g));
end
  
for iter = 1:maxit
  if est
    dir = -H*grad;
    Dp = grad' * dir;
    lambda = [1 0]';
    lambdamin = 2*eps*max(abs(dir)./max(abs(betaold), 1));
    while 1
      if lambda(1) < lambdamin
	beta = betaold;
	break
      end
      beta = betaold + lambda(1)*dir;
      eta = reshape(U*beta, n, g);
      mu = exp(eta - repmat(max(eta, [], 2), 1, g));
      mu = mu./repmat(sum(mu, 2), 1, g);
      if any(any(~mu & G))
	D = inf;
      else
	logmu = G;
	logmu(find(G)) = log(mu(find(G)));
	D = sum(w' * (G .* (logG - logmu) - G + mu));
      end
      if D <= Dold + 1e-4*Dp*lambda(1)
	break
      elseif lambda(1) == 1
	lambda = [-Dp/(2*(D - Dold - Dp)); 1];
      else
	ab = [1 -1; -lambda(2) lambda(1)] * diag(1./lambda.^2) * ...
	     ([D; D2] - Dp*lambda - Dold) / diff(lambda);
	lambda(2) = lambda(1);
	if ab(1) == 0
	  if ab(2) == 0
	    break
	  end
	  lambda(1) = -Dp/(2*ab(2));
	else
	  lambda(1) = (-ab(2) + sqrt(ab(2)^2 - 3*ab(1)*Dp))/(3*ab(1));
	end
      end
      
      if ~isreal(lambda)
	lambda(1) = .1*lambda(2);
      else
	lambda(1) = max(min(lambda(1), .5*lambda(2)), ...
			.1* lambda(2));
      end
      D2 = D;
    end
  elseif g == 2
    deta = eeta./(1 + eeta).^2;
    W = sqrt(w) .* deta ./ sqrt(mu .* (1 - mu));
    z = eta + (G(:,2) - mu)./deta;
    beta = (W(:,ones(p+1,1)) .* U) \ (W .* z);
    eta = U*beta;
    eeta = exp(eta);
    mu = eeta./(1 + eeta);
    logmu = [1-mu, mu];
    if any(any(~logmu & G))
      error('Deviance infinite.')
    end
    logmu(find(G)) = log(logmu(find(G)));
    D = sum(w' * (G .* (logG - logmu)));
  else
    W = sparse(1:n, 1:n, sqrt(w)) * mu ./ sqrt(mu);
    z = eta + (G - mu)./mu;
    beta = (sparse(1:g*n, 1:g*n, W) * U) \ reshape((W .* z), n*g, 1);
    eta = reshape(U*beta, n, g);
    mu = exp(eta);
    if any(any(~mu & G))
      error('Deviance infinite.')
    end
    logmu = G;
    logmu(find(G)) = log(mu(find(G)));
    D = sum(w' * (G .* (logG - logmu) - G + mu));
  end

  if trace
    disp(sprintf('Iter: %d; Dev: %g', iter, 2*D))
  end
  
  if Dold - D < 0
    warning('Deviance diverged.')
    beta = betaold;
    D = Dold;
    break
  elseif Dold - D < D*n*eps
    if trace
      disp('Deviance converged.')
    end
    break
  end
  
  if est
    grad1 = grad;
    delta = mu(:,2:end) - G(:,2:end);
    grad = [delta, (delta(:,repmat((1:g-1)', 1, p)) .* ...
		    X(:, repmat(1:p, g-1, 1)))]' * w;
    dir = beta - betaold;
    if max(dir./max(beta, 1)) < 4*eps
      if trace
	disp('Gradient converged.')
      end
      break
    end
    dg = grad - grad1;
    pdg = dir'*dg;
    Hdg = H*dg;
    gHg = dg'*Hdg;
    u = dir/pdg - Hdg/gHg;
    H = H + dir*dir'/pdg - Hdg*Hdg'/gHg + gHg*u*u';
  end
  
  Dold = D;
  betaold = beta;
end

if ~est & g > 2
  coefs = reshape(beta(1:(g-1)*(p+1)), g-1, p+1);
else
  coefs = reshape(beta, g-1, p+1);
end

if est
  coefs(:,2:p+1) = coefs(:,2:p+1) * diag(1./diff(range));
  coefs(:,1) = coefs(:,1) - coefs(:,2:p+1) * range(1,:)';
end

f = class(struct('coefs', coefs), 'logda', h);  

if nargout > 2
  dev = 2*D;
  
  if nargout > 3
    hess = inv(H);
  end
end





