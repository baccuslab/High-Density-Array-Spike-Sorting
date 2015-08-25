
function [noiseEpochs C] = calculateNoiseEpochsAndCovarianceMatrix(X, spikeEpochs, Tf, varargin)
    % Tf is necessary, in case spikeEpochs is empty!
    P.makeNoiseEpochsShorterBy = 20; % Be sure to not include spikes
    P.minNoiseEpochLength = Tf+100;
    P = mysort.util.parseInputs(P, 'calculateNoiseEpochsAndCovarianceMatrix', varargin);

    
    
    noiseEpochs = mysort.epoch.flip(spikeEpochs, size(X,2));        
    noiseEpochs(:,1) = noiseEpochs(:,1) + P.makeNoiseEpochsShorterBy;
    noiseEpochs(:,2) = noiseEpochs(:,2) - P.makeNoiseEpochsShorterBy;
    C = mysort.util.mcNoiseCovariance(X, noiseEpochs, Tf, P.minNoiseEpochLength);