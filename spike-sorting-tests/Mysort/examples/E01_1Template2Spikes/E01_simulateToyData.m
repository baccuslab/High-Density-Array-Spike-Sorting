function [X, T] = E01_simulateToyData(Tf, nC, distance, noise_sd)
% Creates one simple artificial template and two instances of this
% template. The template will have an amplitude of 1. White Gaussian noise
% is added to the data. Spike sorting is then applied to the data.
% The correct template is used and the result visualised.
%
% Input: 
%     Tf, nC   - The template has nC channels (rows) and is Tf samples long.
%     distance - distance of the start stamples of the two instances of the
%                template                 
%     noise_sd - Standard deviation of the Gaussian noise. 

% Calulate the start samples of the two spikes
offset = 20;
position1 = offset + Tf;
position2 = position1 + distance;

% The length of the whole data
L = offset + position2 + 2*Tf;

% This will be the data
X = zeros(nC,L);

% This is a template that looks like a sine
template_single_channel = sin(linspace(0.02,.98,Tf)*3*pi);
% Create the template in single row notation (channels concatenated)
T = repmat(template_single_channel, 1, nC);
% Create the template in multiple row notation (every row represents one
% channel
T_mat = mysort.util.v2m(T,nC);

% Copy the template in the data
X = mysort.util.shiftSubtract(X,-T_mat,position1);
X = mysort.util.shiftSubtract(X,-T_mat,position2);
% Add some white noise
X = X + noise_sd*randn(size(X));

