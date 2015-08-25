benchmarks = {'Easy1', 'Easy2', 'Difficult1', 'Difficult2'};
noiselevels = {'005', '01', '015', '02'};
dpath = ana.botmpaper.E10_quirogaDataPath();
spike_times = {};
spike_class = {};
nSpikesPerNeuron = [];
benchmarkIdx = [];
noiselevelIdx = [];
%%
count = 0;
for b = 1:length(benchmarks)
    for n = 1:length(noiselevels)
        count = count+1;
        fprintf('Starting with %s %s\n', benchmarks{b}, noiselevels{n});
        quirogaFile = sprintf('C_%s_noise%s.mat', benchmarks{b}, noiselevels{n});      
        D = load(fullfile(dpath, quirogaFile));
        spike_times = [spike_times; D.spike_times];
        spike_class = [spike_class; D.spike_class];
        nSpikesPerNeuron(count,:) = [sum(D.spike_class{1}==1) sum(D.spike_class{1}==2) sum(D.spike_class{1}==3)];
        benchmarkIdx(count) = b;
        noiselevelIdx(count) = n;
    end
end
save('QuirogaBenchmarkDetails.mat', 'spike_times', 'spike_class', 'nSpikesPerNeuron', 'noiselevelIdx', 'benchmarkIdx');