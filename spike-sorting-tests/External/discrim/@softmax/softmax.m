function [f, iter, dev, hess] = softmax(X, k, prior, varargin)
%SOFTMAX Multinomial feed-forward neural-network
%   F = SOFTMAX(X, K, PRIOR) returns a SOFTMAX object containing the
%   weights of a feed-forward neural network trained to minimise the
%   multinomial log-likelihood deviance based on the feature matrix X,
%   class indeces in K and the prior probabilities in PRIOR where
%   PRIOR is optional. See the help for SOFTMAX's parent object class
%   CLASSIFIER for information on the input arguments X, K, and
%   PRIOR. Traditional neural networks minimise the sum squared error,
%   whereas this model assumes that the outputs are a Poisson process
%   conditional on their sum and calculates the error as the residual
%   deviance.
%
%   In addition to the fields defined by the CLASSIFIER class, F
%   contains the following field:
%
%   WEIGHTS: a sparse matrix representing the optimised connection
%   weights where rows represent connections from units that feed
%   other units and columns represent connections to units that are
%   fed by other units. Each non-zero value in this matrix represents
%   a weight connecting unit i to unit j where i is the row and j is
%   the column. There are p+1 input units that are not fed by other
%   units (i.e., the first p+1 columns are all zeros). The first unit
%   always represents the bias, while the following p units represent
%   the inputs to the entire network (i.e., from X). In addition there
%   are g-1 output units where g represents the number of different
%   classes in k. Because output probabilities are normalised over the
%   sum of the exponents, it is assumed that the first class recieves
%   all zero weights and is therefore not explicitly represented in
%   the weight matrix. Output units do not feed other units and
%   therefore there are g-1 less rows than columns (assuming that the
%   missing rows are all zero). All other units are referred to as
%   hidden units and both feed and are fed by other units.
%
%   Because of argument structure ambiguity, PRIOR is not optional
%   when using other options. The default can be assigned by giving
%   an empty PRIOR = [].
%
%   SOFTMAX(X, K, PRIOR, NUNITS, SKIP) where NUNITS is a scalar
%   positive integer and SKIP is either 0 or 1 specifies how many
%   hidden units are present in a single hidden layer neural
%   network. The model is fully connected between adjacent layers. If
%   SKIP is 1, input units are additionally connected to output
%   units. SKIP must be specified when there is only a single hidden
%   layer. If NUNITS is 0, SKIP must also be 0.
%
%   SOFTMAX(X, K, PRIOR, NUNITS) where NUNITS is a vector of positive,
%   non-zero integers of length n specifies how many units are present
%   in each of n hidden layers. All adjacent layers are fully
%   connected, however it is an error to specify a SKIP. If skip
%   wieghts are desired, the weight matrix must be given explicitly
%   (see below).
%
%   SOFTMAX(X, K, PRIOR, WEIGHTS, MASK) where WEIGHTS is a matrix
%   similar to F.WEIGHTS described above, uses the connections and
%   starting weights specified in the matrix. MASK is optional. If
%   given MASK is a matrix the same size as weights consisting of all
%   1's and 0's indicating which weights are to be optimised by the
%   training algorithm. This allows the optimisation of weights that
%   are initially 0 as well as the ability to keep some non-zero
%   weights fixed.
%
%   SOFTMAX(X, K, PRIOR, MASK) is equivalent to SOFTMAX(X, K, PRIOR,
%   WEIGHTS, MASK) where WEIGHTS are assigned randomly. If the initial
%   random weights used by the training algorithm are needed, the MASK
%   argument (or the NUNITS plus SKIP arguments) can be used with a
%   value of 0 for MAXITER (see below).
%
%   By default, SOFTMAX uses no hidden units with skip weights which
%   is functionally equivalent to a logistic discriminant analysis
%   (see LOGDA). However, SOFTMAX will be much slower as the
%   algorithm has been generalised for hidden units.
%
%   SOFTMAX(X, K, PRIOR, ..., DECAY) where DECAY is a positive scalar
%   value less than 1 gives the weight decay for the model. The
%   default decay is 0. DECAY forces the estimate of the residual
%   deviance to be penalised by the magnitude of the estimated
%   weights. Typical values range from .01 for a very large DECAY to a
%   moderate value 10e-6. Because SOFTMAX initially normalises the
%   inputs, this value is independent of the range of X. (However,
%   SOFTMAX rescales the returned weights so that rescaling of input
%   values is not necessary when classifying new data.)
%
%   SOFTMAX(X, K, PRIOR, ..., DECAY, MAXITER) where MAXITER is a
%   positive integer aborts the algorithm after that many
%   iterations. The default value is 200. If a value of 0 is given
%   as MAXITER the algorithm terminates before optimising the
%   connection weights. This is useful for returning a random
%   matrix of weights which can be later manipulated before
%   optimisation. However, if MAXITER is 0, a DECAY value must be
%   given to avoid ambiguity in the arguments.
%
%   SOFTMAX(X, K, PRIOR, ..., MAXITER) is otherwise equivalent to
%   supplying a DECAY of 0 (unless MAXITER is also 0---see above).
%
%   SOFTMAX(X, K, OPTS) allows optional arguments to be passed in the
%   fields of the structure OPTS. Fields that are used by SOFTMAX are
%   PRIOR, NUNITS, SKIPFLAG, WEIGHTS, MASK, DECAY, and
%   MAXITER. However, neither NUNITS nor SKIP may be specified with
%   either WEIGHTS or MASK.
%
%   [F, NITER, DEV, HESS] = SOFTMAX(X, k, ...) Additionally returns
%   the number of iterations required by the algorithm before
%   convergence in NITER, the residual deviance for the fit in DEV and
%   the Hessian matrix of the weights in HESS. HESS is a square matrix
%   where each row and column represents a single weight. The weights
%   are ordered according to the vectorised weight matrix
%   F.WEIGHTS(:);
%
%   SOFTMAX(X, G, ...) where G is a p by g matrix of posterior
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
%   SOFTMAX(F) where F is an object of class LOGDA returns the
%   SOFTMAX equivalent of the logistic discriminant analysis.
%
%   See also CLASSIFIER, LDA, QDA, LOGDA.
%
%   Notes:
%   The argument structure can be rather complicated. The program
%   tries to figure out which argument is which heuristically, but
%   it's probably easy to defeat it. Arguments that are passed to
%   SOFTMAX must be in the order described above although they may
%   be entirely omitted allowing defaults to be used instead.
%
%   References:
%   B. D. Ripley (1996) Pattern Classification and Neural
%   Networks. Cambridge.

%   Copyright (c) 1999 Michael Kiefte.

%   $Id: softmax.m,v 1.1 1999/06/04 18:50:50 michael Exp $
%   $Log: softmax.m,v $
%   Revision 1.1  1999/06/04 18:50:50  michael
%   Initial revision
%

if isa(X, 'logda')
  error(nargchk(1, 1, nargin))
  weights = [sparse(X.nvar+1, X.nvar+1) X.coefs'];
  f = class(struct('weights', weights), 'softmax', X.classifier);  
  return
end

error(nargchk(2, 7, nargin))

if nargin > 2 & isstruct(prior)
  % using option structure
  if nargin > 3
    error(sprintf(['Cannot have arguments following option struct:\n' ...
		   '%s'], nargchk(3, 3, 4)))    
  end
  [prior nhid skip weights mask decay maxit] = ...
      parseopt(prior, 'prior', 'nunits', 'skip', 'weights', 'mask', ...
	       'decay', 'maxiter');  
  if (~isempty(nunits) | ~isempty(skip)) & (~isempty(weights) | ...
					    ~isempty(mask))
    error(['May not specify NUNITS or SKIPFLAG with either WEIGHTS' ...
	   ' or MASK.'])
  end
elseif nargin < 3
  prior = [];
end

[n p] = size(X);

if prod(size(k)) ~= length(k)
  % Multinomial incidence matrix or posterior probabilities
  if length(varargin) > 4
    error(sprintf(['Assuming second argument is an incidence matrix' ...
		   ' of multinomial counts\nor posterior probabilities:' ...
		   ' %s'], nargchk(0, 4, 5)))
  end
  
  [h G w] = classifier(X, k);
  g = size(G, 2);
  logG = G;
  logG(find(G)) = log(G(find(G)));
else
  % Vector of class indeces
  [h G] = classifier(X, k, prior);
  nj = h.counts;
  g = length(nj);
  w = (nj./(n*h.prior))';
  w = w(k);
  logG = 0;
end

% Normalise inputs between (0, 1)
range = h.range;
X = (X - repmat(range(1,:), n, 1)) * diag(1./diff(range));

trace = ~strcmp(warning, 'off');

% varargin will be in this order:
weights = [];
mask = [];
nhid = [];
skip = [];
decay = [];
maxit = [];

if length(varargin)
  % all arguments are real doubles
  if ~isempty(varargin{1}) & isa(varargin{1}, 'double') & ...
	isreal(varargin{1})
    if prod(size(varargin{1})) ~= length(varargin{1})
      %specify weights as matrix
      if length(varargin) >= 2 & ...
	    all(size(varargin{2}) == size(varargin{1}))
	%with mask matrix
	if length(varargin) > 4
	  error(sprintf(['Assuming fifth argument is MASK:' ...
			 ' %s'], nargchk(2, 4, 5)))	  
	end
	varargin = [varargin(1:2), repmat({[]}, 1, 2), varargin(3:end)];
      elseif all(nonzeros(varargin{1}) == 1)
	%only mask matrix
	if length(varargin) > 3
	  error(sprintf(['Assuming fourth argument is MASK:' ...
			 ' %s'], nargchk(1, 3, 4)))
	end
	varargin = [{[]}, varargin(1), repmat({[]}, 1, 2), ...
		    varargin(2:end)];	
      else
	%without mask matrix
	if length(varargin) > 3
	  error(sprintf(['Assuming fourth argument is WEIGHTS:' ...
			 ' %s'], nargchk(1, 3, 4)))
	end
	varargin = [varargin(1), repmat({[]}, 1, 3), ...
		    varargin(2:end)];
      end
    elseif length(varargin{1}) > 1
      % specify number of units in each hidden layer
      if length(varargin) > 3
	error(sprintf(['Assuming fourth argument is the number of' ...
		       ' hidden units\nin each hidden layer:' ...
		       ' %s'], nargchk(1, 3, 4)))
      end
      varargin = [repmat({[]}, 1, 2), varargin(1), {[]}, ...
		  varargin(2:end)];
    elseif round(varargin{1}) == varargin{1}
      if length(varargin) >= 2 & isa(varargin{2}, 'double') & ...
	    isreal(varargin{2}) & length(varargin{2}) == 1 & ...
	    (varargin{2} == 1 | varargin{2} == 0)
	% single hidden layer with skip flag
	if length(varargin) > 4
	  error(sprintf(['Assuming fifth argument is SKIPFLAG:\n' ...
			 ' %s'], nargchk(2, 4, 5)))
	end
	varargin = [repmat({[]}, 1, 2), varargin];
      else
	% third argument is maximum number of iterations
	if length(varargin) > 1
	  error(sprintf(['Assuming fourth argument is MAXITER:' ...
			 ' %s'], nargchk(1, 1, 2)))
	end
	varargin = [repmat({[]}, 1, 5), varargin];
      end
    else
      % third argument is decay
      if length(varargin) > 2
	error(sprintf('Assuming fourth argument DECAY: %s', ...
		      nargchk(1, 2, 3)))	  
      end
      varargin = [repmat({[]}, 4, 1), varargin];
    end
  else
    error('Can''t figure out what third argument should be.')
  end
    
  if length(varargin) == 5 & isa(varargin{5}, 'double') & ...
	length(varargin{5}) == 1 & ...
	round(varargin{5}) == varargin{5}
    % maxiter in decay position
      varargin(5:6) = [{[]}, varargin(5)];
  end
  
  if length(varargin) < 6
    varargin{6} = [];
  end
  
  [weights mask nhid skip decay maxit] = deal(varargin{:});
end

if isempty(decay)
  decay = 0;
elseif ~isa(decay, 'double') | ~isreal(decay) | length(decay) ~= 1 | ...
      decay < 0 | decay >= 1 | isnan(decay)
  error('DECAY must be a positive scalar less than 1.')
end

if ~isempty(weights) | ~isempty(mask)
  normw = 1;
  
  if ~isempty(mask)
    if ~isa(mask, 'double') | ~isreal(mask) | ndims(mask) ~= 2 | ...
	  ~any(any(mask)) | ~all(nonzeros(mask) == 1)
      error('Mask must be a 2-d array of 0s and 1s.')
    end
    
    if isempty(weights)
      weights = mask;
      weights(find(mask)) = 1.4 * rand(nnz(mask), 1) - .7;
      normw = 0;
    elseif ndims(weights) ~= 2 | ~all(size(mask) == size(weights))
      error('MASK and WEIGHTS must be same size.')      
    end
  end

  if normw
    if ~isa(weights, 'double') | ~isreal(weights) | ...
	  ndims(weights) ~= 2 | any(any(isnan(weights)))
      error('Weights must be a 2-d real array.')
    end
    
    if isempty(mask)
      mask = weights;
      mask(find(weights)) = 1;
    end

    % rescale weights because of input normalisation
    weights(1, p+2:end) = weights(1, p+2:end) + ...
	range(1,:) * weights(2:p+1, p+2:end);    
    weights(2:p+1, p+2:end) = diag(diff(range)) ...
	* weights(2:p+1, p+2:end);    
  end
    
  if decay
    % insert redundant output unit to balance weights for decay
    % parameter
    m = sum(weights(:, end-g+2:end), 2)/g;
    weights = [weights(:, 1:end-g+1), -m, ...
	       [weights(:, end-g+2:end) - repmat(m, 1, g-1)]];    
    mask = [mask(:, 1:end-g+1), ...
	    sparse(any(mask(:, end-g+2:end), 2)), ...
	    mask(:, end-g+2:end)];
    % sum of weights from each unit to all outputs now sum to 0
  end
  
  weights(max(find(any(weights | mask, 2))) + 1:end, :) = [];
  mask(size(weights, 1)+1:end, :) = [];
else
  if ~isempty(nhid)
    if ~isa(nhid, 'double') | length(nhid) ~= prod(size(nhid)) ...
	  | ~isreal(nhid) | any(round(nhid) ~= nhid | ...
				isinf(nhid) | nhid < 0)
      error('NUNITS units must be a vector of positive, finite integers.')
    elseif length(nhid) > 1
      if any(nhid <= 0)
	error(['Cannot specify layers with no units with more than one' ...
	       ' hidden layer.'])
      end
    elseif nhid == 0
      nhid = [];
    end
  end

  if isempty(skip)
    skip = 0;
  elseif skip
    if length(nhid) > 1
      error(['May not specify skip weights with more than one hidden' ...
	     ' layer.'])      
    elseif nhid == 0
      error('May not specify skip weights with no hidden layer.')      
    end
    skip = 1;
  else
    skip = 0;
  end
    
  noutput = g - (decay == 0);
  nunits = [p nhid noutput];
  nweights = sum((nunits(1:end-1) + 1) .* nunits(2:end)) + ...
      skip*p*noutput;
  idx = cumsum([2 nunits]);
  nlayer = length(nunits);
  nunits = sum(nunits)+1;
  weights = sparse(nunits - noutput, nunits, nweights);
  mask = weights;
  % connect bias to all hidden and output units.
  mask(1, p+2:end) = 1;
  if skip
    % connect input units to all output units.
    mask(2:p+1, end-noutput+1:end) = 1;
  end
  for i = 1:nlayer-1
    % connect adjacent layers
    mask(idx(i):idx(i+1)-1, idx(i+1):idx(i+2)-1) = 1;
  end
  weights(find(mask)) = 1.4 * rand(nweights, 1) - .7;
end

ninput = max(find(~any(mask | weights)));
if any(~any(mask(:,p+2:end) | weights(:,p+2:end)))
  error(sprintf(['Unit %d has no input:\nInput units must have lowest' ...
		 ' indeces.'], ...
		min(find(~any(mask(:,p+2:end) | ...
			      weights(:,p+2:end))))+p+1)) 
elseif ninput < p+1
  error('Not enough input units.')
elseif ninput > p+1
  error('Too many input units.')
end

noutput = diff(size(mask));
if any(~any(mask | weights, 2))
  error(sprintf(['Unit %d has no ouput:\nOutput units must have' ...
		 ' highest indeces.'], ...
		min(find(~any(mask | weights, 2)))))
elseif noutput < g-(decay == 0)
  error('Not enough output untis.')
elseif noutput > g-(decay == 0)
  error('Not enough input units.')
end

if any(any(tril(weights | mask)))
  error('All weights must join a lower to a higher indexed unit.')
end

if isempty(maxit)
  maxit = 200;
elseif ~isa(maxit, 'double') | ~isreal(maxit) | length(maxit) ~= 1 ...
      | round(maxit) ~= maxit | maxit < 0
  error('MAXITER must be a positive integer.')
end

% initial states
[E post grad] = feedprop(X, G, w, weights, mask, decay);
H = eye(nnz(mask)); % inverse of Hessian
oldweights = weights;
oldE = E;
  
for iter = 1:maxit
  dir = -H*grad; % direction vector
  Ep = grad'*dir; % gradient from vector of partial derivitives (grad)
  lambda = [1 0]'; % length and old length
  lambdamin = 2*eps*max(abs(dir)./max(abs(oldweights(find(mask))), 1));
  while 1
    if lambda(1) < lambdamin
      weights = oldweights;
      break      
    end
    % try point along dir
    weights(find(mask)) = oldweights(find(mask)) + lambda(1)*dir;
    E = feedprop(X, G, w, weights, mask, decay);
    if E <= oldE + 1e-4*Ep*lambda(1)
      break % good enough
    elseif lambda(1) == 1
      lambda = [-Ep/(2*(E - oldE - Ep)); 1];
    else
      ab = [1 -1; -lambda(2) lambda(1)] * diag(1./lambda.^2) * ...
	   ([E; E2] - Ep*lambda - oldE) / diff(lambda);
      lambda(2) = lambda(1);
      if ab(1) == 0
	if ab(2) == 0
	  break
	end
	lambda(1) = -Ep/(2*ab(2));
      else
	labmda(1) = (-ab(2) + sqrt(ab(2)^2 - 3*ab(1)*Ep)) / ...
	    (3* ab(1));
      end
    end
    
    if ~isreal(lambda)
      lambda(1) = .1*lambda(2);
    else
      lambda(1) = max(min(lambda(1), .5*lambda(2)), .1*lambda(2));
    end
    E2 = E;
  end
  
  if trace & ~rem(iter, 10)
    disp(sprintf('Iter: %d; Err: %g', iter, E))
  end
  
  if oldE - E < 0 % indicates divergence (this appears to be normal
                  % for large lambda
    warning('Error diverged.')		  
    weights = oldweights;
    E = oldE;
    break
  elseif oldE - E < E*n*eps % indicates convergence
    if trace
      disp('Error converged.')
    end
    break
  end
  
  grad1 = grad;
  [oldE post grad] = feedprop(X, G, w, weights, mask, decay);
  dir = weights(find(mask)) - oldweights(find(mask));
  if max(dir./max(weights(find(mask)), 1)) < 4*eps
    if trace
      % convergence in weights but not error (probably too strict)
      disp('Gradient converged.')
    end
    break
  end
  oldweights = weights;
  dg = grad - grad1;
  pdg = dir'*dg;
  Hdg = H*dg;
  gHg = dg'*Hdg;
  u = dir/pdg - Hdg/gHg;
  H = H + dir*dir'/pdg - Hdg*Hdg'/gHg + gHg*u*u';  
end

if decay
  % get rid of redundant output by subtracting weights from other
  % output weights
  weights = [weights(:, 1:end-g), [weights(:, end-g+2:end) - ...
		    repmat(weights(:, end-g+1), 1, g-1)]];
end

% rescale to actual non-normalised inputs so no one gets confused.
weights(2:p+1, :) = spdiags(1./diff(range)', 0, p, p) ...
    * weights(2:p+1, :);
weights(1, p+2:end) = weights(1, p+2:end) - ...
    range(1,:) * weights(2:p+1, p+2:end);
f = class(struct('weights', weights), 'softmax', h);

if nargout > 2
  dev = 2*(E - decay*sum(weights(find(mask))));
  
  if nargout > 3
    hess = inv(H);
  end
end

function [E, post, grad] = feedprop(X, G, w, weights, mask, decay)
%FEEDPROP Feed-forward and back-propogate.
%   [E, POST, GRAD] = FEEDPROP(X, G, W, WEIGHTS, MASK, DECAY)
%   returns the error E, the posterior probabilities in the n by g
%   matrix POST and the gradient vector of partial derivitives in
%   the NNZ(MASK) length vector GRAD.

[n p] = size(X);
g = size(G, 2);
logG = G; % so we don't run into NaNs
logG(find(G)) = log(G(find(G)));
ninput = p + 1;
noutput = g - (decay == 0);
nunits = size(weights, 2);
nhidden = nunits - noutput - ninput;

% outputs from individual units including network inputs
y = [ones(n, 1), X, zeros(n, nhidden + noutput)];

for i = ninput+1:nunits % hidden and output units
  [idx, j, wgt] = find(weights(:, i));
  nwgt = length(idx);
  if nwgt % it has happened!
    y(:, i) = sum(y(:, idx) * diag(wgt), 2);
    if i <= nunits - noutput % hidden units are logistic
      out = exp(y(:,i));
      y(:, i) = out./(1+out);
      y(isnan(y(:,i)), i) = 1; % if out is very large
    end
  end
end

% calculate normalised posterior probabilities using SOFTMAX critereon
% if no decay, first class is always zero
post = [zeros(n, decay == 0), y(:, end-noutput+1:end)];
post = exp(post - repmat(max(post(:, 2:end), [], 2), 1, g));
post = post ./ repmat(sum(post, 2), 1, g);
if any(any(~post & G))
  E = inf;
else
  logpost = G; %zeros in G are not a problem. zeros in post mean
		  %that weights are going off into infinity.
  logpost(find(G)) = log(post(find(G)));
  E = sum(w' * (G .* (logG - logpost) - G + post)) + ...
      decay*sum(weights(find(mask)).^2);
end
  
if nargout > 1  
  delta = [zeros(n, nhidden), ...
	   post(:, 1 + (decay == 0):end) - G(:, 1 + (decay == 0):end)];
  
  for i = nhidden:-1:1
    [j idx wgt] = find(weights(i+ninput, :));
    nwgt = length(idx);
    delta(:, i) = y(:, i+ninput) .* (1 - y(:, i+ninput)) .* ...
        sum(delta(:, idx-ninput) * spdiags(wgt', 0, nwgt, nwgt), 2);
  end
  
  [i j] = find(mask);
  grad = (y(:, i) .* delta(:, j-ninput))' * w ...
	 - 2*decay*weights(find(mask));
end




