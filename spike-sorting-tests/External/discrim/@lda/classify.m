function [c, post] = classify(f, X, method, ndims)
%LDA/CLASSIFY Categorise new data with linear discriminants.
%   [C, POST] = CLASSIFY(F, X) classifies the rows of the n by p
%   feature matrix X given the LDA object F, where n is the number of
%   observations or rows in X and p is the number of features or
%   variates. The estimated classes are returned in the length n index
%   vector C, while the posterior probabilities for each class are
%   given in the n by g matrix POST, where g is the number of groups
%   classifiable by F. Each row corresponds to a row in X.
%
%   CLASSIFY(F, X, METHOD) where METHOD is one of 'plug-in',
%   'predictive' or 'debiased' uses the corresponding method for
%   classification. The default is 'plug-in'. The predictive method
%   averages over the uncertainty in the estimation of the mean and
%   transformation matrices using a "vague prior", while the debiased
%   method uses an unbiased estimator of the log probability
%   density. Predictive and debiased estimates are not defined for
%   t-estimator LDA objects. For equal priors in F, all three methods
%   including the plug-in rule for t-parameter estimators, produce the
%   same linear discriminants. When specifying METHOD, only the first
%   few disambiguating letters need be given: i.e., 'pl', 'd' and
%   'pr'.
%
%   CLASSIFY(F, X, METHOD, NDIMS) where ndims is the number of
%   canonical variates to use in the class estimation uses the first
%   ndims columns of the Fisher's linear discriminant transformation
%   matrix returned by LDA/CVAR instead of F.SCALE.
%
%   CLASSIFY(F, X, NDIMS) is equivalent to 
%   CLASSIFY(F, X, 'plug-in', NDIMS).
%
%   CLASSIFY(F, X, OPTS) allows optional arguments to be passed in the
%   fields of the structure options. Fields that are used by
%   LDA/CLASSIFY are NDIMS and METHOD.
%
%   See also LDA, CVAR, CROSSVAL.
%
%   References:
%   B. D. Ripley (1996) Pattern Classification and Neural
%   Networks. Cambridge.

%   Copyright (c) 1999 Michael Kiefte.
%   Additionally based on algorithm presented in S-Plus code written
%   by Ripley and Venables.

%   $Log$

error(nargchk(2, 4, nargin))
classify(f.classifier, X);

M = f.means;
[g p] = size(M);
n = size(X, 1);
S = f.scale;
est = f.est;
nj = f.classifier.counts;
prior = f.classifier.prior;

if nargin > 2 & isstruct(method)
  if nargin > 3
    error(sprintf(['Cannot have arguments following option struct:\n' ...
		   '%s'], nargchk(3, 3, 4)))
  end
  [ndims method] = parseopt(method, 'ndims', 'method');
elseif nargin < 4
  ndims = [];
  if nargin < 3
    method = '';
  elseif ~ischar(method)
    ndims = method;
    method = '';
  end
end

if ~isempty(ndims)
  S = cvar(f);
  if ~isa(ndims, 'double') | length(ndims) ~= 1 | round(ndims) ~= ...
	ndims | ndims < 1
    error(['Number of dimensions NDIMS must be a positive, non-zero' ...
	   ' integer.'])
  elseif ndims > size(S, 2)
    error(sprintf(['Number of dimensions NDIMS can be no more than' ...
		   ' the number of\ncanonical variates (%d)'], ...
		  size(S, 2)))
  end
  S = S(:,1:ndims);
else
  S = f.scale;
  ndims = min(g-1, p);
end

if isempty(method)
  if est == 't'
    method = 't';
  else
    method = 1;
  end
elseif ~ischar(method) | length(method) ~= size(method, 2) | ...
      size(method, 1) ~= 1
  error('METHOD must be a string.');
else
  method = find(strncmp(method, {'plug-in', 'debiased', 'predictive'}, ...
			length(method)));
  if isempty(method)
    error(['METHOD must be one of ''plug-in'', ''debiased'', or' ...
	   ' ''predicitive''.'])
  elseif length(method) > 1
    error('METHOD matches more than one option.')
  end
  
  if est == 't'
    if method == 1
      method = 't';
    else
      error(['May only specify METHOD ''plug-in'' for t-estimator' ...
	     ' LDA objects.'])      
    end
  end
end

if method == 3 & est == 0
  S = S*sqrt(sum(nj)/(sum(nj) - g));
elseif method == 2 & est == 1
  S = S*sqrt((sum(nj) - g)/sum(nj));
end

Mm = mean(M);
Xs = (X - repmat(Mm, n, 1))*S;
Ms = (M - repmat(Mm, g, 1))*S;

switch method
 case 1
  L = repmat(0.5*sum(Ms'.^2) - log(prior), n, 1) - Xs*Ms';
 case 2
  m = sum(nj) - g;
  L = (m - ndims - 1)/(2*m) * (repmat(sum(Ms'.^2)/2, n, 1) - Xs*Ms') ...
      - repmat(log(prior) + p./(2*nj), n, 1);
 case 3
  ntr = sum(nj);
  njadj = nj./(nj + 1);
  L = log(1 + (repmat(sum(Xs.^2, 2), 1, g) - ...
	       2*Xs*Ms' + repmat(sum(Ms'.^2), n, 1)) * ...
	  diag(njadj/ntr)) * (ntr - g + 1)/2 - ...
      repmat(log(prior) + log(njadj)*(p/2), n, 1);
 case 't'
  nu = f.nu;
  L = log(1 + 1/nu*(repmat(sum(Xs.^2, 2), 1, g) - 2*Xs*Ms' + ...
		    repmat(sum(Ms'.^2), n, 1)))*(nu+p)/2 - ...
      repmat(log(prior), n, 1);
end

[Lc c] = min(L, [], 2);

if nargout > 1
  Pr = exp(Lc(:,ones(1, g)) - L);
  post = Pr./repmat(sum(Pr, 2), 1, g);
end







