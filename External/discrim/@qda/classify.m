function [c, post] = classify(f, X, method)
%QDA/CLASSIFY Categorise new data with quadratic discriminants.
%   [C, POST] = CLASSIFY(F, X) classifies the rows of the n by p
%   feature matrix X given the QDA object F, where n is the number of
%   observations or rows in X and p is the number of features or
%   variates. The estimated classes are returned in the length n index
%   vector C, while the posterior probabilities for each class are
%   given in the n by g matrix POST, where g is the number of groups
%   classifiable by F. Each row corresponds to a row in X.
%
%   CLASSIFY(f, X, METHOD) where METHOD is one of 'plug-in',
%   'predictive' or 'debiased' uses the corresponding method for
%   classification. The default is 'plug-in'. The predictive method
%   averages over the uncertainty in the estimation of the mean and
%   transformation matrices using a "vague prior", while the debiased
%   method uses and unbiased estimator of the log probability
%   density. Predictive and debiased estimates are not defined for
%   t-estimator QDA objects. When specifying METHOD, only the first
%   few disambiguating letters need be given: i.e., 'pl', 'd' and
%   'pr'.
%
%   CLASSIFY(f, X, OPTS) allows an optional argument to be passed in
%   the fields of the structure OPTS. The only field that is used by
%   QDA/CLASSIFY is METHOD.
%
%   See also QDA, CROSSVAL.
%
%   References:
%   B. D. Ripley (1996) Pattern Classification and Neural
%   Networks. Cambridge.

%   Copyright (c) 1999 Michael Kiefte.
%   Additionally based on algorithm presented in S-Plus code written
%   by Ripley and Venables.

%   $Log$

error(nargchk(2, 3, nargin))
classify(f.classifier, X);

M = f.means;
[g p]=size(M);
n = size(X, 1);
S = f.scale;
ldet = f.ldet;
est = f.est;
nj = f.classifier.counts;
prior = f.classifier.prior;

if nargin > 2 
  if isstruct(method)
    method = parseopt(method, 'method');
  end
else
  method = '';
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
  S = S.*repmat(shiftdim(sqrt(nj./(nj - 1)), -1), [p p 1]);
  ldet = ldet + p*log((nj - 1)./nj);
elseif method == 2 & est == 1
  S = S.*repmat(shiftdim(sqrt((nj - 1)./nj), -1), [p p 1]);
  ldet = ldet + p*log(nj./(nj - 1));
end

D = zeros(n, g);
for i = 1:g
  D(:,i) = sum(((X - M(i(ones(n, 1)), :)) * S(:,:,i)).^2, 2);
end

switch method
 case 1
  L = 0.5 * D + repmat(ldet/2 - log(prior), n, 1);
 case 2
  if any(nj <= g)
    error(sprintf('Class %d too small.', min(find(nj <= g))))
  end
  m = nj - 1;
  B = p*log(m/2) - sum(digamma((repmat(nj, g, 1) - repmat((1:g)', 1, g))/2));
  L = repmat((ldet + B - p./nj)/2 - log(prior), n, 1) + D * ...
      diag((m - p - 1)./m)/2;
 case 3
  L = repmat(p/2 * log(nj + 1) + gammaln((nj - p)/2) - ...
	     gammaln(nj/2) + ldet/2 - log(prior), n, 1) + ...
      log(1 + D*diag(1./(nj + 1)))*diag(nj/2);
 case 't'
  nu = f.nu;
  L = (nu+p)/2 * log(1 + D/nu) + repmat(ldet/2 - log(prior), n, 1);
end

[Lc c] = min(L, [], 2);
if nargout > 1
  Pr = exp(Lc(:,ones(1, g)) - L);
  post = Pr./repmat(sum(Pr, 2), 1, g);
end
    
function g = digamma(z)

g = zeros(size(z));
if any(any(z < 5))
  x = z(z < 5);
  g(z < 5) = digamma(x + 5) - 1./x - 1./(x + 1) - 1./(x + 2) - 1./(x ...
						  + 3) - 1./(x + 4);
end
if any(any(z >= 5))
  x = z(z >= 5).^-2;
  tail = (x .* (-1/12 + (x .* (1/120 + (x .* (-1/252 + (x .* ...
      (1/240 + (x .* (-1/132 + (x .* (691/32760 + (x .* (-1/12 + ...
      (3617 .* x)/8160))))))))))))));
  g(z >= 5) = log(z(z >= 5)) - 1./(2 * z(z >= 5)) + tail;
end


