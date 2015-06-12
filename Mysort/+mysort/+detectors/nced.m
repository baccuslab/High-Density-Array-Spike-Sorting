
function out = nced(X, p)

% OPTIONAL PARAMETERS
%     p.threshFactor          Determines the threshold of spike detection.
%                             Treshold will be the std of the nced signal
%                             times threshFactor.
%     p.maximaWidth
%     p.method                One of the following
%                                 'sumFirst'  - First sum up the signals of 
%                                             all channels and then square the 
%                                             result to compute the power
%                                 'squareFirst' - First square the signal of
%                                             each channel and then sum up the power

if nargin < 2
    p = [];
end
% make sure the all necessary parameters are set
p = checkNcedParameters(X,p);

% get the number of channels and make sure X has the right format
[nC L] = size(X);
if nC>L
    X = X';
end

if strcmp(p.method,'sumFirst')
    collapseSwitch = 1;
else
    collapseSwitch = 0;
end

% filter the data with "multi channel normalized cumulative energy difference"
ncedSignal = mcNced(X, collapseSwitch);

% do some smoothing and then thresholding to get out the spike times 
out = processSpikeDetectionOutputSignal(ncedSignal, p);

