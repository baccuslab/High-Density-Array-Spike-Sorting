function [CLASSPERF CPUTIMEPERF DETPERF L] = loadQuirogaEval()
benchmarks = {'Easy1', 'Easy2', 'Difficult1', 'Difficult2'};
noiselevels = {'005', '01', '015', '02'};
initTf = 81;
initCutLeft = -5;
Tf = 61;
cutLeft = 15;
dpath = ana.botmpaper.E10_quirogaDataPath();
samplesPerSecond = 32000;


%%
CLASSPERF = [];
CPUTIMEPERF = [];
DETPERF = [];
count = 1;
for b = 1:4 %length(benchmarks)
    for n=1:4 %length(noiselevels)
        savePath = fullfile(dpath, benchmarks{b}, noiselevels{n});
        fprintf('Starting with %s %s\n', benchmarks{b}, noiselevels{n});
        resFile = fullfile(savePath, 'final_results_small');
        L = load(resFile, 'tgdf', 'RES', 'detPerf', 'methods', 'datasetNames');
        correctSpikesIdx = 1:2:length(L.datasetNames);
        for di=1:length(correctSpikesIdx)
            d = correctSpikesIdx(di);
            CLASSPERF(count, di, :) = 100*L.RES(d).performance;
        end
    
        for di=1:length(correctSpikesIdx)
            d = correctSpikesIdx(di)+1;
            CPUTIMEPERF(count, di, :) = mean(L.RES(d).cputime,2);
        end              
        DETPERF(count,:,:) = -100*L.detPerf;
        count = count+1;
    end
end
fprintf('Loading done.\n');