% Example 1 a) "the use of the sorter"
% See E01_readme for details.

% Set the template length
Tf = 10;
% Set the number of channels
nC = 2;
% Set the distance of the two spikes
dist = 15;
% Set the noise standard deviation
noise_sd = .08;
% Generate a piece of toy data (X) with two spikes (T is the template)
[X, T] = E01_simulateToyData(Tf, nC, dist, noise_sd);

C = (noise_sd^2)*eye(Tf*nC);
NE = mysort.util.NoiseEstimator(C, Tf);

% Init the sorter
botm = mysort.sorters.BOTM(NE, Tf, T, 'upsample', 10);

% Sort our artifical data
gdf = botm.sort(X);

% And plot the piece of data
botm.plotLastChunkSorting('figureTitle', sprintf('Example1 dist=%d', dist));

