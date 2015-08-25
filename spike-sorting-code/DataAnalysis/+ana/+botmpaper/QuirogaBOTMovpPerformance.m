benchmarks = {'Easy1', 'Easy2', 'Difficult1', 'Difficult2'};
noiselevels = {'005', '01', '015', '02'};
initTf = 81;
initCutLeft = -5;
Tf = 61;
cutLeft = 15;
dpath = ana.botmpaper.E10_quirogaDataPath();
samplesPerSecond = 32000;
Eucl  = @ana.botmpaper.Eucl;
Maha  = @ana.botmpaper.Maha;
Conv  = @ana.botmpaper.Conv;
Match = @ana.botmpaper.Match;
Botm  = @ana.botmpaper.Botm;
global savefig__;
savefig__ = 1;
figSavePath = 'C:\Dropbox\PaperBOTM1\DetectorFigures\';
cols = {'b', 'g', 'r'};
%%
RR = [];
b = 1 %length(benchmarks)
n = 1 %length(noiselevels)
savePath = fullfile(dpath, benchmarks{b}, noiselevels{n});
if ~exist(savePath, 'file')
    mkdir(savePath);
end      

pic = 0;
fprintf('Starting with %s %s\n', benchmarks{b}, noiselevels{n});
close all
quirogaFile = sprintf('C_%s_noise%s.mat', benchmarks{b}, noiselevels{n});        
GT = ana.botmpaper.E10_preprocessing(dpath, quirogaFile, initCutLeft, initTf, dpath);

bufferFile = fullfile(savePath, 'botmOutputDistributionsThrehold.mat');
   
% RUN QuirogaOptimalThreshold first to save all buffered files.
load(bufferFile);
%% COMPUTE FILTER OUTPUTS
Y = zeros(size(GT.X));
Tf2 = floor(size(T,2)/2);
M = zeros(2*Tf2+1, 3, 3);
for i=1:3
    t = T(i,:);
    f = t/C_;
    Y(i,:)  = conv(GT.X, flipud(f(:)), 'same') - .5*t*f';
    for j=1:3
        M(:,i,j) = conv(T(j,:), flipud(f(:)), 'same');
    end
end
maxTau = 8;
DC = mysort.sorters.DiscriminantFunctionContainer(Y', M, maxTau);
% d1 = DC.get(1,2,1);
% d2 = DC.get(2,1,-1);
% figure
% plot(d1)
% hold on
% plot(d2,'g')
%%
figure
plot(T', 'linewidth', 3)
figure
plot(Y(:,1900:2200)', 'linewidth', 3)
hold on
plot(DC.DOvp(1900:2200,:))

%%
[maxY maxSingleTemp] = max(Y);

maxYOvp      = zeros(3, size(Y,2));
maxYOvpClass = zeros(3, size(Y,2));
peaksOvp = {};
peakIds = [];
for i=1:3
    myidx = DC.DOvpIndex(:,2)==i;
    myIds = DC.DOvpIndex(myidx,1);
    myShifts = DC.DOvpIndex(myidx,4);
    myBoarderIds = myIds(abs(myShifts)==maxTau);
    
    % First, let the overlapping functions compete with each other
    Ytmp = DC.DOvp(:,myidx)';
    [maxYOvp(i,:) maxYOvpClass(i,:)] = max(Ytmp);
    
    % now find the local maxima
    [pksO, locsO] = findpeaks(maxYOvp(i,:), 'MINPEAKDISTANCE', ceil(Tf/2), 'MINPEAKHEIGHT', log(.99));
    
    % if the local maxima was an overlapping function at -maxTau or
    % +maxTau, we dont know if it was actually tau or bigger than tau. so
    % we ignore it
    ids = maxYOvpClass(i,locs0);
    pks0(ismember(ids, myBoarderIds)) = [];
    locs0(ismember(ids, myBoarderIds)) = [];
%     Y(i,:); 
    
    peaksOvp{i} = [locsO(:) pksO(:) myIds(maxYOvpClass(i,locsO)')];
    peakIds(i,:) = myIds';
end

%%
allPeaks = sortrows(cell2mat(peaksOvp'));
keepIdx = zeros(size(allPeaks,1),1);
% find spikes
count = 1;
while count <=size(allPeaks,1)
    k = count+1;
    while (k<=size(allPeaks,1)) && ( (allPeaks(k,1) - allPeaks(count,1)) < 20)
        k=k+1;
    end
    [mx mxidx] = max(allPeaks(count:k-1,2));
    keepIdx(count+mxidx-1) = 1;
    count=k;
end
allPeaks = allPeaks(keepIdx==1,:);

%
[gdf_ovp wasResolved] = mysort.spiketrain.resolveOvpIndexInGdf(allPeaks(:,[3 1]), DC.DOvpIndex);
gdf_ovp(:,2) = gdf_ovp(:,2)-cutLeft;

R_ovp = mysort.spiketrain.alignGDFs(tgdf, gdf_ovp, 10, 3, 25);
mysort.plot.printEvaluationTable(R_ovp)
%%
def = mysort.util.defs();
fpidx = find(R_ovp.spikeLabel2{1} == def.fp);
fp_times = R_ovp.St2{1}(fpidx);
for k=1:5
    fpt = fp_times(k);
    xrange = fpt-400:fpt+400;
    apidx      = ((fpt-400) < allPeaks(:,1)) & (allPeaks(:,1) < (fpt+400));
    tgdfidx    =  (fpt-400) < tgdf(:,2)   &    tgdf(:,2) < (fpt+400);
    gdf_ovpidx =  (fpt-400) < gdf_ovp(:,2)   &    gdf_ovp(:,2) < (fpt+400);
    
    figure
    plot(xrange,Y(:,xrange)')
    hold on
    for i=1:3
        idx = ((fpt-400) < peaksOvp{i}(:,1)) & (peaksOvp{i}(:,1) < (fpt+400));
        plot(peaksOvp{i}(idx,1), peaksOvp{i}(idx,2), 'o', 'markersize', 14, 'color', cols{i});
    %     mysort.plot.vectorColor(peaksOvp{i}(:,3)))
    end

    plot(allPeaks(apidx,1), allPeaks(apidx,2), 'cx', 'markersize', 12, 'linewidth', 2)
    
    mysort.plot.spiketrain(tgdf(tgdfidx,:), 'yoffset', -100, 'colorList', cols)
    mysort.plot.spiketrain(gdf_ovp(gdf_ovpidx,:), 'colorList', cols)
    plot(allPeaks(apidx,1), allPeaks(apidx,2), 'cx', 'markersize', 12, 'linewidth', 2)
    mysort.plot.dataCursorMode();
    axis tight
end

%%
figure
plot(Y')
hold on
for i=1:3
    plot(peaksOvp{i}(:,1), peaksOvp{i}(:,2), 'o', 'markersize', 14, 'color', cols{i});
%     mysort.plot.vectorColor(peaksOvp{i}(:,3)))
end
plot(allPeaks(:,1), allPeaks(:,2), 'cx', 'markersize', 12, 'linewidth', 2)
mysort.plot.dataCursorMode();
%% 
ovpgdf = [];
for i=1:3
    for k=1:size(peaksOvp{i},1)
        mySpikeTime   = peaksOvp{i}(k,1);
        mySpikeHeight = peaksOvp{i}(k,2);
        mySpikeId     = peakIds(i,peaksOvp{i}(k,3));
        keepMe = true;
        for i2=1:3
            if i2==i
                continue
            end
            otherPeaks = find(abs(peaksOvp{i2}(:,1) - mySpikeTime) < Tf/2);
            if isempty(otherPeaks)
                % no other spike here
            elseif peaksOvp{i2}(otherPeaks,1) < mySpikeHeight
                % other spike is smaller
            else
                keepMe = false;
            end 
        end
        if keepMe
            ovpgdf(end+1,:) = [mySpikeId mySpikeTime];
        end
    end
end
[gdf_ovp wasResolved] = mysort.spiketrain.resolveOvpIndexInGdf(ovpgdf, DC.DOvpIndex);
gdf_ovp(:,2) = gdf_ovp(:,2)-cutLeft;
% %%
% gdf_ovp = gdf_single;
% counts = [0 0 0];
% for i=1:size(peaksO_clean,1)
%     idx = find(abs(gdf_single(:,2)-(peaksO_clean(i,1)-cutLeft)) < 5);
%     ovp = [maxOvpTemp(peaksO_clean(i,1))+3 peaksO_clean(i,1)-cutLeft];
%     if isempty(idx)
%         gdf_ovp(end+1, :) = ovp;
%         counts(1) = counts(1)+1;
%     elseif peaksO_clean(i,2) > peaksY(idx,2)
%         gdf_ovp(idx, :) = ovp;
%         counts(2) = counts(2)+1;
%     else
%         counts(3) = counts(3)+1;
%     end
% end
% 
% counts
% [gdf_ovp wasResolved] = mysort.spiketrain.resolveOvpIndexInGdf(gdf_ovp, DC.DOvpIndex);
R_ovp = mysort.spiketrain.alignGDFs(tgdf, gdf_ovp, 10, 3, 25);
mysort.plot.printEvaluationTable(R_ovp)

% [tgdf(1:10,:) gdf_single(1:10,:)]
% figure
% plot(Ytmp')
% look where the overlap functions at the border are maximal. those we
% ignore since it could be, that the overlap is actually further apart and
% we dont want to detect the same spike twice.
%% Find peaks
[pksY, locsY] = findpeaks(maxY, 'MINPEAKDISTANCE', 10, 'MINPEAKHEIGHT', log(.99));
peaksY = [locsY(:) pksY(:)];

%% Block overlap functions at the boarders
% borderOvpIdx = ismember(DC.DOvpIndex(:,4), [-maxTau maxTau]);
% borderOvpIDs = DC.DOvpIndex(borderOvpIdx,1);
% peaksO_clean = peaksO;
% peaksO_clean(ismember(maxOvpTemp(locsO), borderOvpIDs-3),:) = [];
% size(peaksO,1) - size(peaksO_clean,1)
% size(peaksO_clean,1)

%% GDF single spikes
gdf_single = [maxSingleTemp(peaksY(:,1))' peaksY(:,1)-cutLeft];
R_single = mysort.spiketrain.alignGDFs(tgdf, gdf_single, 5, 3, 25);
mysort.plot.printEvaluationTable(R_single)
% [tgdf(1:10,:) gdf_single(1:10,:)]

%% GDF overlapping spikes

figure
plot(gdf_single(:,2), zeros(1, size(gdf_single,1)), 'xr', 'markersize', 14, 'linewidth', 3)
hold on
plot(peaksO_clean(:,1)-cutLeft, ones(1, size(peaksO_clean,1)), 'xb',  'markersize', 14, 'linewidth', 3)
set(gca, 'ylim', [-10 10], 'xlim', [100 5000]);



%%
figure;
plot(maxY)
hold on
plot(maxOvp, 'r')


%% annotation of peaks
OVERLAP_DIST = 15;
OTHERSPIKE_DIST = 6;
SINGLESPIKE_DIST = 3;

