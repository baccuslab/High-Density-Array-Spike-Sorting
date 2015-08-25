close all
clear all
PREC = ana.botmpaper.harrisAnalysis.preprocess();
MS   = ana.botmpaper.harrisAnalysis.meanshiftSort();
PRECGT = ana.botmpaper.harrisAnalysis.preprocessGroundTruth();
CLUEVAL = ana.botmpaper.harrisAnalysis.clusterEval();
VN = PREC.getValidTrialNames();
global savefig__;
savefig__ = 0;
%%
% R = PREC.run('d11221.002');
% PREC.makeFigures('d11221.002')
% % 

%%
PREC.resetTrial('d533102');
MS.resetTrial('d533102');

%%
PREC.resetTrial('d11221.002');
MS.resetTrial('d11221.002');

%%
MS.resetTrial('d11222.001');

%%
cd(PREC.figurePath)
cd(PREC.dataOutPutPath)

RMS = MS.run('d11221.002');
MS.makeFigures('d11221.002')

%%
MS.runAllValidTrials()

%%
PRECGT.resetTrial('d14521.002');
PRECGT.run('d11221.002');
PRECGT.makeFigures('d14521.002')
%%
PRECGT.makeAllFigures()
%%
PRECGT.runAllValidTrials()

%%
CLUEVAL.resetTrial('d11221.002');
CLUEVAL.run('d11221.002');
%%
CLUEVAL.resetTrial('d11222.001');
CLUEVAL.run('d11222.001');
%%
CLUEVAL.resetTrial('d12821.001');
CLUEVAL.run('d12821.001');
%%
CLUEVAL.makeFigures('d11221.002')
%%
trialName = 'd11221.002';
PREC.resetTrial(trialName);
PRECGT.resetTrial(trialName);
MS.resetTrial(trialName);
CLUEVAL.resetTrial(trialName);

PREC.run(trialName);
PRECGT.run(trialName);
MS.run(trialName);
CLUEVAL.run(trialName);

PREC.makeFigures(trialName);
PRECGT.makeFigures(trialName);
MS.makeFigures(trialName);
CLUEVAL.makeFigures(trialName);

%%
CLUEVAL.runAllValidTrials()
CLUEVAL.makeAllFigures()
% CLUEVAL.resetAllInProcess()
%-CLUEVAL.resetAllTrials()
