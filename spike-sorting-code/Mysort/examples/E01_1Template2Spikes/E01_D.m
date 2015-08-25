% Example 1 d) "signal to noise ratio and multiple channels"
% See E01_readme for details.


% This is basically the same example as E01_C but we increase the number of
% channels and not the sample rate.
Tf = 10;
nC = 1;
noise_sd = 1;
dist = 100;
[X, T] = E01_simulateToyData(Tf_1, nC, dist, noise_sd);


C = (noise_sd^2)*eye(Tf*nC);
NE = mysort.util.NoiseEstimator(C, Tf);
% Init the sorter
botm = mysort.sorters.BOTM(NE, Tf, T);

botm.sort(X);    
botm.plotLastChunkSorting('figureTitle', sprintf('Example 1 D) nC=%d', nC));

% Note, that no spike is visible and nothing is detected.


% Now we create the same data, however we increase the number of channels.
% The signal to noise ratio (if defined by the peak of the template and the
% noise standard deviation alone) will stay the same.
nC = 14;
[X, T] = E01_simulateToyData(Tf_1, nC, dist, noise_sd);

C = (noise_sd^2)*eye(Tf*nC);
NE = mysort.util.NoiseEstimator(C, Tf);
% Init the sorter
botm = mysort.sorters.BOTM(NE, Tf, T);

botm.sort(X);    
botm.plotLastChunkSorting('figureTitle', sprintf('Example 1 D) nC=%d', nC));

% Note, that the two spikes can be clearly detected, highlighting the fact,
% that the "peak to noise sd ratio" is not a good definition for the signal to
% noise ratio. This effect is of course only so dramatic, if the noise is
% not strongly correlated between channels.
