%% INIT
pd = pdefs();
expPath = fullfile(pd.serverData, '..', 'Mea1k', 'shared', '140821', 'data', 'hamster_DS_4');

mapFile = fullfile(pd.serverData, '..', 'Mea1k', 'shared', '140821', 'map.mat');
rawFile = fullfile(expPath, '0000.raw.h5');

outFolder = 'C:\LocalData\tmp';
if ~exist(outFolder, 'file'); mkdir(outFolder); end
M = load(mapFile);

deflationr        = [0 1 2];
chunkSizer        = [0 200 1000 10000];
chunkEachChannelr = [0 1];

combis = mysort.util.buildCombinatorics([length(deflationr) length(chunkSizer) length(chunkEachChannelr)]);

%% TEST SUITE
for i = 1:size(combis,1)
    deflation = deflationr(combis(i,1));
    chunkSize = chunkSizer(combis(i,2));
    chunkEachChannel = chunkEachChannelr(combis(i,3));
    
    convertedFile = fullfile(outFolder, sprintf('TestFile_chunkEachChannel_%d_chunkSize_%d_deflation_%d.h5', chunkEachChannel, chunkSize, deflation));
    if exist(convertedFile, 'file');
        continue;
    end
%     clear DS; delete(convertedFile); end

    mysort.mea.preprocessMea1kH5File(rawFile, M.map,...
                            'prefilter', 1, 'outFile', convertedFile, ...
                            'subtractMeanOverAllChannels', 1, ...
                            'restrictToTimePeriod', [10000 80000], ...
                            'chunkLen', chunkSize, ...
                            'deflation', deflation,...
                            'chunkIndividualChannels', chunkEachChannel);
    
end

%% SPECIAL CASES
convertedFile = fullfile(outFolder, sprintf('TestFile_noChunks.h5'));
if ~exist(convertedFile, 'file');
    mysort.mea.preprocessMea1kH5File(rawFile, M.map,...
                        'prefilter', 1, 'outFile', convertedFile, ...
                        'subtractMeanOverAllChannels', 1, ...
                        'restrictToTimePeriod', [10000 80000], ...
                        'chunkLen', 0);
end

