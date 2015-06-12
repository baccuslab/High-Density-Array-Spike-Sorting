function e = shrink(f, gamma)
%LDA/SHRINK Shrink within-groups covariance matrix to identity.
%   E = SHRINK(F, GAMMA) shrinks the within-groups covariance matrix
%   to a diagonal matrix proportional to EYE(p) where p is the number
%   of features in the trianing set of F. GAMMA must be a real
%   positive value no more than 1. A GAMMA of 1 corresponds to
%   complete shrinkage while a GAMMA of 0 corresponds to no shrinkage.
%
%   See also LDA, COV.
%
%   References: 
%   B. D. Ripley (1996) Pattern Classification and Neural
%   Networks. Cambridge.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

error(nargchk(2, 2, nargin))

if isempty(gamma) | ~isa(gamma, 'double') | ~isreal(gamma) | ...
      length(gamma) ~= 1 | gamma < 0 | gamma > 1 | isnan(gamma)
  error(['Scale parameter GAMMA must be a positive scalar no more' ...
	 ' than 1.'])
end

e = f;
if gamma
  s = inv(f.scale);
  c = s'*s;
  p = size(c, 1);
  e.scale = inv(chol((1 - gamma)*c + gamma/p*trace(c)*eye(p)));  
end




