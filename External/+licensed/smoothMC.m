

function smoothX = smoothMC(X, span, method)
% Smoothes the rows in X with given span and method

if nargin < 2
    span = 3;
end
if nargin < 3
    method = 'sgolay';
end



[nrSamples sampleLength] = size(X);
smoothX = zeros(nrSamples, sampleLength);

for i=1:nrSamples
   smoothX(i,:) = smooth(X(i,:), span, method);
end
    