function [c, post, y] = classify(f, X, opt)
%SOFTMAX/CLASSIFY Categorise new data with SOFTMAX object.
%   [C, POST] = CLASSIFY(F, X) classifies the rows of the n by p
%   feature matrix X given the LOGDA object F, where n is the number
%   of observations or rows in X and p is the number of features or
%   variates. The estimated classes are returned in the length n index
%   vector C, while the posterior probabilities for each class are
%   given in the n by g matrix POST with each row corresponding to a
%   row in X.
%
%   [C, POST, OUT] = CLASSIFY(F, X) additionally returns the
%   outputs for individual units in the network in matrix
%   POST. POST has as many rows as X with a column for each unit in
%   the network including the bias unit and the original input
%   matrix X (i.e., the first p+1 columns are [ZEROS(n,1) X]). The
%   remaining columns represent the outputs of the hidden and
%   output layers.
%
%   See also SOFTMAX, CROSSVAL.
%
%   References:
%   B. D. Ripley (1996) Pattern Classification and Neural
%   Networks. Cambridge.

%   $Id: classify.m,v 1.1 1999/06/04 18:50:50 michael Exp $
%   $Log: classify.m,v $
%   Revision 1.1  1999/06/04 18:50:50  michael
%   Initial revision
%

%   Copyright (C) 1999 Michael Kiefte.

error(nargchk(2, 3, nargin))
classify(f.classifier, X);

weights = f.weights;
g = diff(size(weights))+1;
p = sum(~any(weights))-1;
n = size(X, 1);

ninput = p+1;
noutput = g-1;
nunits = size(weights, 2);
nhidden = nunits - noutput - ninput;

y = [ones(n, 1), X, zeros(n, nhidden + noutput)];

for i = ninput+1:nunits
  [idx, j, wgt] = find(weights(:, i));
  nwgt = length(idx);
  if nwgt
    y(:, i) = sum(y(:, idx) * diag(wgt), 2);
    if i <= nunits - noutput
      out = exp(y(:,i));
      y(:, i) = out./(1+out);
      y(isnan(y(:,i)), i) = 1;
    end
  end
end

post = [zeros(n, 1), y(:, end-noutput+1:end)];
post = exp(post - repmat(max(post(:, 2:end), [], 2), 1, g));
post = post ./ repmat(sum(post, 2), 1, g);
[m c] = max(post, [], 2);








