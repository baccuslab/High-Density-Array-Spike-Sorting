% Example 1 c) "signal to noise ratio and the template length"
% See E01_readme for details.


% In this example we create first a piece of data, that is so heavily
% corrupted by noise, that we cannot detect the spikes anymore:
Tf_1 = 10;
Tf_2 = 100;
nC = 1;
noise_sd = 1;
dist = 2*Tf_2;

[X, T] = E01_simulateToyData(Tf_1, nC, dist, noise_sd);

% Init the sorter
C = (noise_sd^2)*eye(Tf*nC);
NE = mysort.util.NoiseEstimator(C, Tf);

% Init the sorter
botm = mysort.sorters.BOTM(NE, Tf, T, 'upsample', 10);

% Sort our artifical data
botm.sort(X);    
fprintf('Example 1 b) Distance = %d\n', dist);

% And plot the piece of data
botm.plotLastChunkSorting('figureTitle', sprintf('Example1 C) Tf=%d', Tf_1));

% Note, that no spike is visible and nothing is detected.


% Now we create the same data, however we increase the number of channels.
% The signal to noise ratio (if defined by the peak of the template and the
% noise standard deviation alone) will stay the same.
C = (noise_sd^2)*eye(Tf_2*nC);
NE = mysort.util.NoiseEstimator(C, Tf_2);
[X, T] = E01_simulateToyData(Tf_2, nC, dist, noise_sd);

botm = mysort.sorters.BOTM(NE, Tf_2, T, 'upsample', 10);

botm.sort(X);    
botm.plotLastChunkSorting('figureTitle', sprintf('Example1 C) Tf=%d', Tf_2));

% Note, that the template can be clearly detected, highlighting the fact,
% that the "peak to noise sd ratio" is not a good definition for the signal to
% noise ratio. This effect is of course only so dramatic, if the noise is
% white and Gaussian.
