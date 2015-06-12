
function [sorting SO] = sort(X, varargin)
    % myosort.sort(X) performs spike sorting on the data in X.
    % INPUT:
    %  X   - either a matrix containing multichannel data (rows correspond
    %        single recording channels) or a dataObject
    %  smap_per_sec  (optional) - samples per second for the recording in X
    %
    % OUTPUT:
    %  sorting - a matrix with two columns. Every row corresponds to one 
    %            identified spike in the recordings in X. The first column
    %            indicates the neuron ID, the second the time point for the
    %            spike
    %  SO  - sorting object. A class containing detailed information about
    %        the spike sorting and functions to access those details and
    %        visualize them
    P = mysort.configuration();
    
    P.method = 'MeanSpike';
    P.Tf = 55;             % template length in samples
    P.maxElNoiseDist = 50; % microns, maximal distance of two electrodes
                           % for which the noise covariance needs to be computed
    P.maxNoiseLag = 15;    % samples, maximal lag for which the noise cov
                           % needs to be computed
    P.thr_noise = 3.5;     % below is noise in x*std
    P.thr_spikes = 4;      % above are spikes in x*std
    P = util.parseInputs(P, 'mysort.sort', varargin);
    
    % Input checks
    X = datasource.DataObject(X, varargin{:});
    
    % Calculate signal noise standard deviation
    [SO.smad, SO.epochs_noise] = noise.estimateSigma(X, 2*P.Tf, P.thr_noise);
    SO.epochs_spikes = noise.detectSpikeEpochs(X, P.thr_spikes*SO.smad, P.Tf);
    
    % Estimate noise covariance matrix
    SO.C = noise.Covest(X, 'noiseEpochs', noiseEpochs, ...
        'maxDist', P.maxElNoiseDist, 'maxLag', P.maxNoiseLag);
    SO.R = noise.Covest(X, 'maxDist', P.maxElNoiseDist, 'maxLag', P.maxNoiseLag);
    
    % Switch for method
    possibleMethods = {'MeanSpike'};
    
    for i=1:length(possibleMethods)
        if strcmp(P.method, possibleMethods{i})
            eval(possibleMethods{i})
            return
        end
    end
    possibleMethods = [possibleMethods repmat(' - ', size(possibleMethods,1), 1)]';
    error('Method: %s not known! Choose one of %s', P.method, [possibleMethods{:}]);
    
    
    function MeanSpike()
        disp('bla');
    end
end
    
    