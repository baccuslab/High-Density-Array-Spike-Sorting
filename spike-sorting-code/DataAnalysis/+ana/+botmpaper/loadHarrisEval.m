function [REJECTIONPERF, CLASSPERF, DETPERF, CPUTIMEPERF, correctSpikesIdx, otherSpikesIdx, L, nTemplatesPerDataSet] = loadHarrisEval()

    REJECTIONPERF = [];
    CLASSPERF = [];
    DETPERF = [];
    CPUTIMEPERF = [];
    %% Load Data            
    names = [1:4 6];
    %% DO NOT USE THE 10kHz FILE !!
    names = [1:4];
    nTemplatesPerDataSet = [];
    for nameidx = 1:length(names)
     name = names(nameidx);
     fprintf('Starting with %d\n', name);
        fprintf('Loading %d.\n', name);
        [D.path D.outpath] = ana.harris.datapath();
        D.info = ana.harris.info(name);
        D.folder   = D.info{1};
        D.name     = D.info{2};
        D.srate    = D.info{3};
        D.nC       = length(D.info{7});
        D.Tf       = D.info{8};
        D.cutLeft  = D.info{9};
        D.nUnits   = D.info{11};   
        dpath = fullfile(D.path, '..', 'MeanShiftSorting6', D.name);

        
        L = load(fullfile(dpath, 'final_results_small'));
        nTemplatesPerDataSet(nameidx) = size(L.RES(1).TM,3);
        nD = length(L.datasetNames);
        correctSpikesIdx = 2:3:nD;
        otherSpikesIdx   = 1:3:nD-1;        
        for di=1:length(correctSpikesIdx)
            d = correctSpikesIdx(di);
            for m=1:size(L.methods,1)
                CLASSPERF(nameidx, di, m) = 100*L.RES(d).classifications(L.correctTemplateIdx,m)/sum(L.RES(d).classifications(:,m));
            end
        end
        for di=1:length(otherSpikesIdx)
            d = otherSpikesIdx(di);
            for m=1:size(L.methods,1)
                REJECTIONPERF(nameidx, di, m) = (100-100*L.RES(d).classifications(L.correctTemplateIdx,m)/sum(L.RES(d).classifications(:,m)));
            end
        end        
        for di=1:length(otherSpikesIdx)
            d = otherSpikesIdx(di);
            for m=1:size(L.methods,1)
                CPUTIMEPERF(nameidx, di, m) = mean(L.RES(d).cputime(m,:));
            end
        end              
        DETPERF(nameidx,:,:) = -100*L.detPerf;
    end
    fprintf('Loading done.\n');