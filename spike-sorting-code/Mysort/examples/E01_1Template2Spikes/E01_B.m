% Example 1 b) "bursts, dead time, overlap with myself"
% See E01_readme for details.

% Create a simple example in which the two spike get closer and closer.
Tf = 10;
nC = 1;
noise_sd = .1;
for dist=[10 8 5 4]
    [X, T] = E01_simulateToyData(Tf, nC, dist, noise_sd);

    C = (noise_sd^2)*eye(Tf*nC);
    NE = mysort.util.NoiseEstimator(C, Tf);

    % Init the sorter
    botm = mysort.sorters.BOTM(NE, Tf, T, 'upsample', 1);
    
    % Sort our artifical data
    botm.sort(X);    
    fprintf('Example 1 b) Distance = %d\n', dist);
    % And plot the piece of data
    botm.plotLastChunkSorting('figureTitle', sprintf('Example1 dist=%d', dist));
end

% For this template we can see, that even of its length of 10 samples, we
% can still resolve the overlap of 5 samples. For 4 sample only one spike
% might be detected, depending on the noise.