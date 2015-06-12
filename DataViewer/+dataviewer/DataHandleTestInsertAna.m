%%
DB = dataviewer.DataHandle('dbconfig', db.munk.DBconfig_INSERT, ...
    'log_function', @disp);

%%
monkey = 'Julia';
exp = '5143';
block = 'a';
tetrodeNr = 1;
trialIdxList = 1:5;
gdf = [ 1 100
        2 200
        3 300
        4 400
        5 500];
algoNeuronIDs = 1:5;
gdfList = {gdf, gdf, gdf, gdf, gdf};

T = randn(90, 4, 5);
cutleft = 20;
C = xcorr(randn(1000,4), 89);

algoName = 'test';

% DB.DBH.queryFieldnames('analysis')
% DB.query('SELECT * FROM analysis LIMIT 20');

%%
DB.insertSorting(monkey, exp, block, tetrodeNr, trialIdxList, algoName, algoNeuronIDs, gdfList, T, cutleft, C);

tId = 18227;
anaID = 1435;
CC = DB.getCovarianceMatrix('trialIDs', tId, 'analysisIDs', anaID);
CC = mysort.noise.xcov2ccol(CC);
figure;
subplot(1,2,1);
imagesc(mysort.noise.xcorr2Cte(C));
subplot(1,2,2);
imagesc(mysort.noise.ccol2Cte(CC));
