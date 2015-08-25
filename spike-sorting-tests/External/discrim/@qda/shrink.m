function e = shrink(f, alpha, gamma)
%QDA/SHRINK Covariance shrinkage.
%   E = SHRINK(F, ALPHA) shrinks the covariances matrix of each class
%   in the QDA object F toward the pooled covariance matrix as
%   returned by COV. ALPHA must be a real positive value no more than
%   1. An ALPHA of 1 corresponds to complete shrinkage and effectively
%   converts F to a linear discriminant analysis while an ALPHA of 0
%   corresponds to no shrinkage.
%
%   E = SHRINK(F, ALPHA, GAMMA) first shrinks the covariance
%   matrices toward the pooled covariance matrix and then shrinks
%   the results toward a diagonal matrix proportional to the
%   identity. A GAMMA of 1 corresponds to complete shrinkage while
%   a GAMMA of 0 corresponds to none.
%
%   References: 
%   B. D. Ripley (1996) Pattern Classification and Neural
%   Networks. Cambridge.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

error(nargchk(2, 3, nargin))

if isempty(alpha) | ~isa(alpha, 'double') | ~isreal(alpha) | ...
      length(alpha) ~= 1 | alpha < 0 | alpha > 1 | isnan(alpha)
  error(['Scale parameter ALPHA must be a positive scalar no'  ...
	 ' greater than 1.'])
end

if nargin < 3
  gamma = 0;
elseif isempty(gamma) | ~isa(gamma, 'double') | ~isreal(gamma) | ...
      length(gamma) ~= 1 | gamma < 0 | gamma > 1 | isnan(gamma)
  error(['Scale parameter GAMMA must be a positive scalar no greater' ...
	 ' than 1.'])  
end
  
e = f;
if alpha | gamma
  nj = f.classifier.counts;
  n = sum(nj);
  [g p] = size(f.means);
  
  C = zeros(p, p, g);
  for i = 1:g
    S = inv(f.scale(:, :, i));
    C(:,:,i) = S'*S;
  end

  if alpha
    S = sum(repmat(shiftdim(nj/n, -1), [p p 1]) .* C, 3);
    for i = 1:g
      C(:,:,i) = ((1 - alpha)*nj(i)*C(:,:,i) + alpha*n*S) / ...
	  ((1 - alpha)*nj(i) + alpha*n);
    end
  end
  if gamma
    trace = sum(reshape(C(repmat((0:p-1)'*p + (1:p)', 1, g) + ...
			  repmat(p*p*(0:g-1), p, 1)), p, g));
    C(:, :) = (1 - gamma)*C(:, :) + sparse(repmat((1:p)', g, 1), ...
					   1:p*g, reshape(repmat((gamma/p) ...
						  * trace', p, 1), p*g, 1));
  end
  for i = 1:g
    r = chol(C(:,:,i));
    e.ldet(i) = 2*sum(log(abs(diag(r))));
    e.scale(:, :, i) = inv(r);
  end
end





