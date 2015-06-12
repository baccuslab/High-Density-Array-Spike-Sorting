function [lambda, ratio] = cvar(f)
%LDA/CVAR Fisher's linear discriminant analysis.
%   [LAMBDA, RATIO] = CVAR(F) return the canonical variates of the
%   Fisher's linear discriminant analysis in the columns of LAMBDA
%   based on the LDA object F. LAMBDA is a scaling matrix which
%   maximises the ratio of the between to within-groups
%   variance. RATIO gives the percentage of between-group variance
%   accounted for by the corresponding canonical variate.
%
%   References:
%   B. D. Ripley (1996) Pattern Classification and Neural
%   Networks. Cambridge.   

%   Copyright (c) 1999 Michael Kiefte.
%   Based also on algorithm presented in S-Plus code written by Ripley
%   and Venables.   

%   $Log$

prior = f.classifier.prior;
M = f.means;
S = f.scale;
g = length(prior);
n = sum(f.classifier.counts);
K = f.est ~= 1;
Ms = diag(sqrt(n*prior/(g - K))) * ...
     (M - repmat(prior * M, g, 1)) * S;
[u s v] = svd(Ms, 0);
r = sum(diag(s) > n*eps*s(1));
lambda = S * v(:, 1:r);

if nargout > 1
  s = diag(s(1:r, 1:r)).^2;
  ratio = s/sum(s);
end
