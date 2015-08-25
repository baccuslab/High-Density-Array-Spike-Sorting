% Testscript for the botm initialization procedure. Simulates a simple test
% scenario and computed the filter coefficients in a semi-efficient way -
% on a subset of all electrodes for each neuron, depending on the neurons
% template. 
%% INIT
L = 200000;
nC = 100;

nSpikes = 200;
nT = 10;
T = [0 .5 .75 .8 .75 .3 -1 -15 -10 -5 -3 -1.5 0 1 1.5 1.75 1.8 1.75 1.5 1 .75 .5 .3 .15 .1 .05 0];
T = T/norm(T);

p_gt_zero = .1;

spikeTrains = cell(1, nT);
connectionMatrix = zeros(nT, nC);

%% Compute Noise
X = 5*randn(L,nC);
% correlate "neighboring" electrodes
for i=1:nC-1
    X(:,i) = (X(:,i) + .1*X(:,i+1))/1.1;
end
X(:,1) = (X(:,1) + .1*X(:,end))/1.1;
% correlate first and 6th channel
X(:,1) = (X(:,1) + .5*X(:,6))/1.5;
X(:,6) = (X(:,6) + .5*X(:,1))/1.5;

%% Build connection matrix between electrodes and neurons and copy a spike into the noise 
for i=1:nT
    spikeTrains{i} = round(rand(1, nSpikes)*(L-200) + 100);
    connectedElectrodes = find(rand(1, nC)<p_gt_zero);
    amplitudes = 20*randn(1, length(connectedElectrodes)).^2;
    connectionMatrix(i,connectedElectrodes) = amplitudes;
    for j=1:length(spikeTrains{i})
        for k=1:length(amplitudes)
            tidx = (spikeTrains{i}(j):spikeTrains{i}(j)+length(T)-1)-10;
            X(tidx,connectedElectrodes(k)) = X(tidx,connectedElectrodes(k)) + ...
                                             amplitudes(k)*T';
        end
    end
end

%% Smooth noise and spikes
X = filter([1 1 1 1 1 1 1]/7, 1, X);
% Wrap the data in a datasource object
DS = mysort.ds.Matrix(X, 20000);

%% ------------------------------------------------------------------------
% END OF SETTING UP THE SIMULATION - DO THE INIT COMPUTATION
% -------------------------------------------------------------------------
% DS  - A data source object that must give access to the extracellular data
%       with the same preprocessing filter as will be used by the hardware.
% spiketrains  - The spike trains for each neuron (a cell array, one cell
%                element per neuron) containing the spike trains that were
%                established by the initial spike sorting.
% COMPUTE THE ACTUAL INITIALIZATION
BIO = mysort.mea.BOTMInitOptimizer(DS, spikeTrains, ...
                    'maxNSpikesForTemplateComputation', 500,...
                    'blockSpikesForNoiseComputationThreshold', 1);
% And save to the init file for the HW Spike Sorter interface
BIO.saveHWSpikeSorterInitFile('BOTMInitOptimizerTestSave.h5');


%% ------------------------------------------------------------------------
% THE REST IS TO CHECK IF THE COMPUTATIONS WERE OK
% -------------------------------------------------------------------------

% Compute the botm discriminants
spike_prior = 0.01;
Y = zeros(size(X,1), nT);
CONF = zeros(BIO.P.filterLength,nT);
for k=1:nT
    nCeffective = BIO.connectionMatrix(k,1);
    myChannels  = BIO.connectionMatrix(k,2:nCeffective+1) + 1;
    F = mysort.wf.v2t(BIO.filterCoefficients(k,:), nCeffective);
    T = mysort.wf.v2t(BIO.effectiveTemplates(k,:), nCeffective);
    for j=1:nCeffective
        idx = myChannels(j);
        Y(:,k)    = Y(:,k)    + conv2(X(:,idx), flipud(F(:,j)), 'same');
        CONF(:,k) = CONF(:,k) + conv2(T(:,j), flipud(F(:,j)), 'same');
    end
    Y(:,k) = Y(:,k) - .5*BIO.botmConstants(k) + log(spike_prior);
end

figure, ah = subplot(2,1,1); plot(X(1:10000,:)), title('Raw data'), ah(2) = subplot(2,1,2); plot(Y(1:10000,:)), linkaxes(ah, 'x');
figure, plot(CONF), title('Responses of Filter on Targets');

%% VISUALIZATION
C = BIO.XCorrContainer.getCte4Channels(1:10);

figure, imagesc(C), title('Covariance Matrix (for some electrodes)');
figure, plot(BIO.filterCoefficients'), title('Filter Coefficients (on selected electrodes)');
figure, plot(BIO.effectiveTemplates'), title('Effective Templates (on selected electrodes)');
figure, plot(mysort.wf.t2v(BIO.templates)'), title('Templates (on all electrodes)');

