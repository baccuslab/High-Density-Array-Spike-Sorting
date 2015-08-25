
function [data trialStart trialEnd] = readATF(fileName)
% reads .atf files
try
    dataStartRow = 11;
    delimeter = '\t';

    nrChannels = 6;

    s = importdata(fileName, delimeter, dataStartRow);
    data = s.data;
    trialLength = size(data,1);
    nrTrials = (size(data,2)-1) / nrChannels;
    expLength = trialLength * nrTrials;

    trialStart = 1:trialLength:expLength;
    trialEnd = trialLength:trialLength:expLength;

    % TODO: what info does the first column represent???
    % remove the first column of the data matrix
    v = data(:,1);
    data(:,1) = [];

    % reshape the data such that the rows are channels and the colums are time
    % points
    channelIdxStart = 1 : nrChannels : nrTrials*nrChannels;
    channelIdxEnd = nrChannels : nrChannels : nrTrials*nrChannels;
    nData = zeros(nrChannels, expLength);
    for i=1:nrTrials
        nData(:, trialStart(i):trialEnd(i)) = data(:, channelIdxStart(i):channelIdxEnd(i))';
    end

    data = nData;
catch ME
    fprintf('Error while reading file %s', fileName);
    rethrow(ME)
end
