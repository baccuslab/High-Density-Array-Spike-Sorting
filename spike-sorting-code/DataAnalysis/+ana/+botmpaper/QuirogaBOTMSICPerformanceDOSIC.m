function [gdf Y1] = QuirogaBOTMSICPerformanceDOSIC(Y, M, refactory, refactoryPeaks, NOISEPRIOR, SPIKEPRIOR)
nT = size(M,2);
Tf2 = floor(size(M,1)/2);

%% FIRST DETECTION
[Dmax MaxIds] = max(Y, [], 1);
spikeEpochs = mysort.epoch.fromBinaryVectorMinLen(...
                        Dmax>log(NOISEPRIOR), 2*refactory+1);
L = mysort.epoch.length(spikeEpochs);

locs = zeros(length(L),1);
for i=1:length(L)
    [mx mxidx] = max(Dmax(spikeEpochs(i,1):spikeEpochs(i,2)));
    locs(i) = spikeEpochs(i,1) + mxidx - 1;
end

% % make long epochs short
% spikeEpochs(L>2*refactory+1,2) = spikeEpochs(L>2*refactory+1,1) + 2*refactory;
% nE = size(spikeEpochs,1);
% % [pks, locs] = findpeaks(Dmax, 'MINPEAKDISTANCE', refactoryPeaks, 'MINPEAKHEIGHT', log(NOISEPRIOR));
% 
% %% MAKE SURE PEAKS ARE REALLY MAXIMA INSIDE WINDOW, FINDPEAKS IS NOT OPTIMAL IN THIS RESPECT
% idx = repmat(0:2*refactory, nE, 1) + repmat(spikeEpochs(:,1), 1, 2*refactory+1);
% idx = idx';
% idx(idx<1) = 1;
% idx(idx>size(Dmax,2)) = size(Dmax,2);
% 
% X = reshape(Dmax(idx(:)), 2*refactory+1, nE);
% [mx mxidx] = max(X,[], 1);
% locs = spikeEpochs(:,1) + mxidx(:) - 1;


%% MAKE GDF
gdf = [MaxIds(locs(:))' locs(:)];

%% SIC
Y1 = Y;

for i=1:nT
    sub = squeeze(M(:,:,i))' + log(SPIKEPRIOR);
    sub(i,:) = inf;
    mySpikeTimes = gdf(gdf(:,1)==i,2);
    mySpikeTimes(mySpikeTimes <= Tf2) = [];
    mySpikeTimes(mySpikeTimes >= size(Y1,2)-Tf2) = [];
    if ~isempty(mySpikeTimes)
        IDX = repmat(mySpikeTimes, 1, size(M,1)) + repmat( (-Tf2:Tf2), length(mySpikeTimes), 1);

        for k=1:length(mySpikeTimes)
            Y1(:,IDX(k,:)) = Y1(:,IDX(k,:)) - sub;
        end
    end
end

%     else
%         k = 1;
%         offset = 300;
%         figure;
%         ah = subplot(2,1,1);
%         plot(Y1(:,  mySpikeTimes(k)-offset:mySpikeTimes(k)+offset)', 'linewidth', 2)
%         hold on
%         xrange = offset+1-Tf2:offset+Tf2+1;
%         plot(xrange, sub');
%         
%         Y1 = mysort.util.shiftSubtract(Y1, sub, mySpikeTimes(k)-Tf2-1);
%         ah(2) =subplot(2,1,2);
%         plot(Y1(:,  mySpikeTimes(k)-offset:mySpikeTimes(k)+offset)', ':')
%         linkaxes(ah, 'xy');
%     end