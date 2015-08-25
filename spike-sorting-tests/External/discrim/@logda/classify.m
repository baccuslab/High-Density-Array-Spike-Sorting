function [c, post] = classify(f, X, opt)
%LOGDA/CLASSIFY Categorise new data with logistic discriminants.
%   [C, POST] = CLASSIFY(F, X) classifies the rows of the n by p
%   feature matrix X given the LOGDA object F, where n is the number
%   of observations or rows in X and p is the number of features or
%   variates. The estimated classes are returned in the length n index
%   vector C, while the posterior probabilities for each class are
%   given in the n by g matrix POST, where g is the number of groups
%   classifiable by F. Each row corresponds to a row in X.
%
%   See also LOGDA, CROSSVAL.
%
%   References:
%   P McCullagh and J. A. Nelder (1989) Generalized Linear
%   Models. Second Edition. Chapman & Hall.
%   B. D. Ripley (1996) Pattern Classification and Neural
%   Networks. Cambridge.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

error(nargchk(2, 3, nargin))
classify(f.classifier, X);

beta = f.coefs;
g = size(beta, 1) + 1;
p = size(beta, 2) - 1;
n = size(X, 1);

col = sparse(n+1:g*n, repmat(1:g-1, n, 1), 1, g*n, g-1);
U = [col, [col(:,repmat((1:g-1)', 1, p)) .* ...
	   repmat(X(:,repmat(1:p, g-1, 1)), g, 1)]];
L = reshape(U*beta(:), n, g);
[Lc c] = max(L, [], 2);

if nargout > 1
  Pr = exp(L - Lc(:,ones(1, g)));
  post = Pr./repmat(sum(Pr, 2), 1, g);
end









