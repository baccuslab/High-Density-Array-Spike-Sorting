%% sortingDemonstration
% Choose an existing data path for tempory data output
close all
% clear all
tmpDataPath = 'C:\Users\frankef\Desktop\sortingTest12\';
% tmpDataPath = '/home/frankef/bel.svn/hima_internal/cmosmea/trunk/matlab/testsorting';

%% Make some simple toy data
nSamplesPerSecond = 20000;
nChannel = 4;
nSeconds = 10;
N = randn(nSeconds*nSamplesPerSecond, nChannel);

% the noise is only white so far. Lets add some covariance between the
% channels:
varPerChannel = [30 35 40 32];
noiseCorrelationMatrix = [1  .2 .1 .1
                          .2 1  .2 .1
                          .1 .2 1  .2
                          .1 .1 .2 1];

C = zeros(nChannel, nChannel);
for c1 = 1:nChannel
    C(c1, c1) = varPerChannel(c1);
    for c2 = c1+1:nChannel
        C(c1, c2) = noiseCorrelationMatrix(c1, c2)*sqrt(varPerChannel(c1)*varPerChannel(c2));
        C(c2, c1) = C(c1, c2);
    end
end
N = N*chol(C);
% Now we give the noise also some correlation in time and a somewhat more
% realistic frequency composition
hpf = 300;  %Hz
lpf = 6000; %Hz
forder = 101;
b = mysort.mea.filter_design_fir(hpf, lpf, nSamplesPerSecond, forder);
N = conv2(N, b', 'same');

% So far we have noise but no spikes. So lets add some spikes.
t1 = -[1:6 5:-1:-6 -5:1:0]; t1 = 25*t1/abs(max(t1));
T1 = [    [t1   0  0    0]
          [0 .4*t1 0    0]
          [0    0 .2*t1 0]
          [0    0  0    .1*t1]]';

t2 = -cos(linspace(0, pi, length(t1)-4)); t2 = 50*t2/abs(max(t2));
T2 = [    [.3*t2 0    0  0]
          [0 .2*t2    0  0]
          [0     0 .5*t2 0]
          [0     0    0 .7*t2]]';

margin = 100; %samples
firingRate1 = 30; %Hz
firingRate2 = 40; %Hz
spikeTimes1 = sort(randi(nSeconds*nSamplesPerSecond-2*margin, nSeconds*firingRate1, 1)+margin);
spikeTimes2 = sort(randi(nSeconds*nSamplesPerSecond-2*margin, nSeconds*firingRate2, 1)+margin);
% remove refractory period violations
spikeTimes1([false; diff(spikeTimes1)<size(T1,1)]) = [];
spikeTimes2([false; diff(spikeTimes2)<size(T2,1)]) = [];
spikeTimes1([false; diff(spikeTimes1)<size(T1,1)]) = [];
spikeTimes2([false; diff(spikeTimes2)<size(T2,1)]) = [];

groundTruth = sortrows([[ones(length(spikeTimes1), 1); 2*ones(length(spikeTimes2), 1)] [spikeTimes1; spikeTimes2]], 2);
idx1 = (repmat(spikeTimes1,1,size(T1,1)) + repmat(1:size(T1,1), length(spikeTimes1), 1))';
idx2 = (repmat(spikeTimes2,1,size(T2,1)) + repmat(1:size(T2,1), length(spikeTimes2), 1))';

% Create some amplitude variation
A1 = 1 + .1*randn(length(spikeTimes1), 1);
A2 = 1 + .1*randn(length(spikeTimes2), 1);

% Now create our simulated spikes
S = zeros(size(N));
S(idx1(:),:) = S(idx1(:),:) + reshape(permute(mysort.wf.v2t(A1*mysort.wf.m2v(T1'), nChannel), [1 3 2]), numel(idx1), nChannel);
S(idx2(:),:) = S(idx2(:),:) + reshape(permute(mysort.wf.v2t(A2*mysort.wf.m2v(T2'), nChannel), [1 3 2]), numel(idx2), nChannel);
S = conv2(S, b', 'same');

% And our data as a superposition of spikes and noise
X = (N+S);

%% if you want to cluster only a single channel, uncomment next line
% X = X(:, 1); 

%% PLOT DATA with Ground truth spike trains superimposed
DS_noiseOnly = mysort.ds.Matrix(N, nSamplesPerSecond, 'Noise only');
DS_spikesOnly = mysort.ds.Matrix(S, nSamplesPerSecond, 'Spikes only');
DS = mysort.ds.Matrix(X, nSamplesPerSecond, 'ToyData');
SSC = mysort.spiketrain.SpikeSortingContainer('GroundTruth', groundTruth, 'wfDataSource', DS);
DS.addSpikeSorting(SSC);
mysort.plot.SliderDataAxes({DS_noiseOnly, DS_spikesOnly, DS}, 'channelSpacer', repmat(3*sqrt(mean(varPerChannel)), 3, 1));

%% SORT DATA
%--------------------------------------------------------------------------
% Set the dead-time of the spike detector to the length of the spikes
P = struct();
P.spikeDetection.minDist = size(T1,1);   
% Set the time in which multiple detection on different channels are merged
% into one spike to the length of the templates
P.spikeDetection.mergeSpikesMaxDist = size(T1,1); 
% This controls degree of clustering. Set to 1.1 to get more clusters, set to 1.4 to get less clusters
P.clustering.meanShiftBandWidthFactor = 1.1;
P.botm.run = 1;                          % activate botm sorting
%--------------------------------------------------------------------------
%---  This is the call to the actual sorting function:  -------------------
[S P] = mysort.sorters.sort(DS, tmpDataPath, 'testSorting', P);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
clusterSorting = [S.clusteringMerged.ids(:) S.clusteringMatched.ts(:)];
botmSorting    = S.botm.gdf;

%% Plot Sorting In Cluster Space, Cut and Aligned Spikes
mysort.plot.waveformsVertical(S.spikeCut.wfs, 'IDs', S.clusteringMerged.ids, 'stacked', 0, 'forceV', 1, 'nC', size(DS,2));
mysort.plot.figureTitle('Toy Data, Cut Waveforms');

mysort.plot.waveformsVertical(S.spikeAligned.wfs, 'IDs', S.clusteringMerged.ids, 'stacked', 0, 'forceV', 1, 'nC', size(DS,2));
mysort.plot.figureTitle('Toy Data, Aligned Waveforms');

mysort.plot.clustering(S.spikeFeatures.X, S.clusteringMerged.ids);
mysort.plot.figureTitle('Toy Data, Cluster Space (PCA)');

%% Plot Sorting On Continuous Data
DS_clustering = mysort.ds.Matrix(X, nSamplesPerSecond, 'Cluster-Sorting');
SSC2 = mysort.spiketrain.SpikeSortingContainer('ClusterSorting', clusterSorting, 'wfDataSource', DS);
DS_clustering.addSpikeSorting(SSC2);

DS_botm = mysort.ds.Matrix(X, nSamplesPerSecond, 'BOTM-Sorting');
SSC3 = mysort.spiketrain.SpikeSortingContainer('BOTMSorting', botmSorting, 'wfDataSource', DS);
DS_botm.addSpikeSorting(SSC3);

mysort.plot.SliderDataAxes({DS, DS_clustering, DS_botm}, 'channelSpacer', repmat(3*sqrt(mean(varPerChannel)), 3, 1));

%% Compute Performance
maxJitter = 10; % samples
maxShift = 20;  % samples
maxOverlapDist = 20; % samples
Rclu = mysort.spiketrain.alignGDFs(groundTruth, clusterSorting, maxJitter, maxShift, maxOverlapDist);
Rbotm = mysort.spiketrain.alignGDFs(groundTruth, S.botm.gdf, maxJitter, maxShift, maxOverlapDist);
mysort.plot.printEvaluationTable(Rclu);
mysort.plot.printEvaluationTable(Rbotm);
fprintf('Total number of spike sorting errors\n')
fprintf(' Clustering: %d\n', sum(Rclu.totErr));
fprintf(' BOTM      : %d\n', sum(Rbotm.totErr));

%% Analyze errors
% mysort.plot.sortingErrors(Rbotm, S.botmsorter, 'FNO');