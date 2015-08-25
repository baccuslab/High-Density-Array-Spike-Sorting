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
runName = 'V01';
outPath = '/links/groups/hima/recordings/collaborations/Quiroga/BOTMPaper1OvpRun/';
outPathRun = fullfile(outPath, runName);
if ~exist(outPath, 'file'); mkdir(outPath); end
if ~exist(outPathRun, 'file'); mkdir(outPathRun); end
RR = [];
lambda = .3;
resampleP = 3;
resultT = zeros(4*4, 3);
i = 0;
for b = 1:length(benchmarks)
    for n = 1:length(noiselevels)
        i = i+1;
        fname = fullfile(outPathRun, sprintf('OvpRunR_%d_%d.mat', b, n));
        if ~exist(fname, 'file')
            ana.botmpaper.QuirogaBOTMovpPerformance2Body()
            save(fname, 'R_ovp', 'gdf_resolved', 'tgdfshifted', 'gdf_raw', 'gdf_ovp', 'Y', 'DMax', 'Ids', 'lambda', 'resampleP')
        else
            disp('Loading...')
            load(fname);
        end
        resultT(i,:) = [sum(R_ovp.detErr) sum(R_ovp.nCL) sum(R_ovp.totErr)];
    end
end
resultT
return

%% COMPILE RESULTS
resultT = zeros(4*4, 3);
i = 0;
for b = 1:length(benchmarks)
    for n = 1:length(noiselevels)
        i=i+1;
        fname = fullfile(outPathRun, sprintf('OvpRunR_%d_%d.mat', b, n));
        try
            D = load(fname, 'R_ovp', 'gdf_resolved', 'tgdfshifted', 'gdf_raw', 'gdf_ovp', 'Y', 'DMax', 'Ids');
            resultT(i,:) = [sum(D.R_ovp.detErr) sum(D.R_ovp.nCL) sum(D.R_ovp.totErr)];
        end
    end
end    
resultT

return
%%
b = 1
n = 1
fname = fullfile(outPathRun, sprintf('OvpRunR_%d_%d.mat', b, n));
load(fname, 'R_ovp', 'gdf_resolved', 'tgdfshifted', 'gdf_raw', 'gdf_ovp', 'Y', 'DMax', 'Ids');
mysort.plot.printEvaluationTable(R_ovp)

%% PLOT FPs OVERLAPS !
def = mysort.util.defs();
stidx = 1;
fpidx = find(R_ovp.spikeLabel2{stidx} == def.fp);
fp_times = R_ovp.St2{stidx}(fpidx);
for k=1:min(8,length(fp_times))
    fpt = fp_times(k);
    xrange = max(1,fpt-400):fpt+400;
    tgdfidx    =  (fpt-400) < tgdfshifted(:,2)   &    tgdfshifted(:,2) < (fpt+400);
    gdfidx =  (fpt-400) < gdf_ovp(:,2)   &    gdf_ovp(:,2) < (fpt+400);
    gdf_rawidx = (fpt-400) < gdf_raw(:,2)   &    gdf_raw(:,2) < (fpt+400);
    
    figure
    ah = subplot(3,1,1);
    plot(xrange,Y(:,xrange)')
    hold on
%     for i=1:3
%         idx = ((fpt-400) < peaksOvp{i}(:,1)) & (peaksOvp{i}(:,1) < (fpt+400));
%         plot(peaksOvp{i}(idx,1), peaksOvp{i}(idx,2), 'o', 'markersize', 14, 'color', cols{i});
%     %     mysort.plot.vectorColor(peaksOvp{i}(:,3)))
%     end

   
    mysort.plot.spiketrain(tgdfshifted(tgdfidx,:), 'yoffset', -400, 'colorList', cols)
    mysort.plot.spiketrain(gdf_ovp(gdfidx,:), 'colorList', cols)
%     plot(allPeaks(apidx,1), allPeaks(apidx,2), 'cx', 'markersize', 12, 'linewidth', 2)
    mysort.plot.dataCursorMode();
    axis tight
    
    ah(2) = subplot(3,1,2);
    plot(xrange, DMax(xrange,:))
    hold on
    plot(gdf_raw(gdf_rawidx,2), DMax(gdf_raw(gdf_rawidx,2),:), 'x', 'markersize', 13', 'linewidth',2);
    % plot ovps
    ovps_time = gdf_raw(gdf_rawidx,2);
    ovps_ =  Ids(gdf_raw(gdf_rawidx,2),:)>3;
    [mx mxid] = max(DMax(gdf_raw(gdf_rawidx,2),:),[],2);
    ovps = [];
    for kk=1:length(mxid)
        ovps(kk) = ovps_(kk,mxid(kk));
    end
    plot(ovps_time(ovps==1), mx(ovps==1), 'oc', 'markersize', 13', 'linewidth',2);
    
    ah(3) = subplot(3,1,3);
    plot(xrange, Ids(xrange,:));
    hold on
    plot(gdf_raw(gdf_rawidx,2), Ids(gdf_raw(gdf_rawidx,2),:), 'x', 'markersize', 13', 'linewidth',2);
        
    linkaxes(ah, 'x'); 
    mysort.plot.dataCursorMode();
end

%% PLOT FNs OVERLAPS !
def = mysort.util.defs();
stidx = 3;
fpidx = find(R_ovp.spikeLabel1{stidx} == def.fn);
fp_times = R_ovp.St1{stidx}(fpidx);
for k=1:min(3,length(fp_times))
    fpt = fp_times(k);
    xrange = fpt-400:fpt+400;
    tgdfidx    =  (fpt-400) < tgdfshifted(:,2)   &    tgdfshifted(:,2) < (fpt+400);
    gdfidx =  (fpt-400) < gdf_ovp(:,2)   &    gdf_ovp(:,2) < (fpt+400);
    gdf_rawidx = (fpt-400) < gdf_raw(:,2)   &    gdf_raw(:,2) < (fpt+400);
    
    figure
    ah = subplot(3,1,1);
    plot(xrange,Y(:,xrange)')
    hold on
%     for i=1:3
%         idx = ((fpt-400) < peaksOvp{i}(:,1)) & (peaksOvp{i}(:,1) < (fpt+400));
%         plot(peaksOvp{i}(idx,1), peaksOvp{i}(idx,2), 'o', 'markersize', 14, 'color', cols{i});
%     %     mysort.plot.vectorColor(peaksOvp{i}(:,3)))
%     end

   
    mysort.plot.spiketrain(tgdfshifted(tgdfidx,:), 'yoffset', -400, 'colorList', cols)
    mysort.plot.spiketrain(gdf_ovp(gdfidx,:), 'colorList', cols)
%     plot(allPeaks(apidx,1), allPeaks(apidx,2), 'cx', 'markersize', 12, 'linewidth', 2)
    mysort.plot.dataCursorMode();
    axis tight
    
    ah(2) = subplot(3,1,2);
    plot(xrange, DMax(xrange,:))
    hold on
    plot(gdf_raw(gdf_rawidx,2), DMax(gdf_raw(gdf_rawidx,2),:), 'x', 'markersize', 13', 'linewidth',2);
    % plot ovps
    ovps_time = gdf_raw(gdf_rawidx,2);
    ovps_ =  Ids(gdf_raw(gdf_rawidx,2),:)>3;
    [mx mxid] = max(DMax(gdf_raw(gdf_rawidx,2),:),[],2);
    ovps = [];
    for kk=1:length(mxid)
        ovps(kk) = ovps_(kk,mxid(kk));
    end
    plot(ovps_time(ovps==1), mx(ovps==1), 'oc', 'markersize', 13', 'linewidth',2);
    
    ah(3) = subplot(3,1,3);
    plot(xrange, Ids(xrange,:));
    hold on
    plot(gdf_raw(gdf_rawidx,2), Ids(gdf_raw(gdf_rawidx,2),:), 'x', 'markersize', 13', 'linewidth',2);
        
    linkaxes(ah, 'x'); 
    mysort.plot.dataCursorMode();
end


%% PLOT CLASSIFICATION ERRORs OVERLAPS !
def = mysort.util.defs();
stidx = 3;
fpidx = find(R_ovp.spikeLabel1{stidx} == def.cl);
fp_times = R_ovp.St1{stidx}(fpidx);
for k=1:min(3,length(fp_times))
    fpt = fp_times(k);
    xrange = fpt-400:fpt+400;
    tgdfidx    =  (fpt-400) < tgdfshifted(:,2)   &    tgdfshifted(:,2) < (fpt+400);
    gdfidx =  (fpt-400) < gdf_ovp(:,2)   &    gdf_ovp(:,2) < (fpt+400);
    gdf_rawidx = (fpt-400) < gdf_raw(:,2)   &    gdf_raw(:,2) < (fpt+400);
    
    figure
    ah = subplot(3,1,1);
    plot(xrange,Y(:,xrange)')
    hold on
%     for i=1:3
%         idx = ((fpt-400) < peaksOvp{i}(:,1)) & (peaksOvp{i}(:,1) < (fpt+400));
%         plot(peaksOvp{i}(idx,1), peaksOvp{i}(idx,2), 'o', 'markersize', 14, 'color', cols{i});
%     %     mysort.plot.vectorColor(peaksOvp{i}(:,3)))
%     end

   
    mysort.plot.spiketrain(tgdfshifted(tgdfidx,:), 'yoffset', -400, 'colorList', cols)
    mysort.plot.spiketrain(gdf_ovp(gdfidx,:), 'colorList', cols)
%     plot(allPeaks(apidx,1), allPeaks(apidx,2), 'cx', 'markersize', 12, 'linewidth', 2)
    mysort.plot.dataCursorMode();
    axis tight
    
    ah(2) = subplot(3,1,2);
    plot(xrange, DMax(xrange,:))
    hold on
    plot(gdf_raw(gdf_rawidx,2), DMax(gdf_raw(gdf_rawidx,2),:), 'x', 'markersize', 13', 'linewidth',2);
    % plot ovps
    ovps_time = gdf_raw(gdf_rawidx,2);
    ovps_ =  Ids(gdf_raw(gdf_rawidx,2),:)>3;
    [mx mxid] = max(DMax(gdf_raw(gdf_rawidx,2),:),[],2);
    ovps = [];
    for kk=1:length(mxid)
        ovps(kk) = ovps_(kk,mxid(kk));
    end
    plot(ovps_time(ovps==1), mx(ovps==1), 'oc', 'markersize', 13', 'linewidth',2);
    
    ah(3) = subplot(3,1,3);
    plot(xrange, Ids(xrange,:));
    hold on
    plot(gdf_raw(gdf_rawidx,2), Ids(gdf_raw(gdf_rawidx,2),:), 'x', 'markersize', 13', 'linewidth',2);
        
    linkaxes(ah, 'x'); 
    mysort.plot.dataCursorMode();

end

return

%%
figure
ah = subplot(2,1,1);
plot(DMax)
ah(2) = subplot(2,1,2);
plot(Ids);
linkaxes(ah, 'x');

%% Find peaks
pks = {};
for i=1:3
    [pksY, locsY] = findpeaks(Y(i,:), 'MINPEAKDISTANCE', ceil(Tf/2), 'MINPEAKHEIGHT', log(.99));
    pks{i,1} = [locsY(:)-cutLeft pksY(:) ones(length(pksY), 1)*i];
end
allPeaks = sortrows(cell2mat(pks));
%%
keepIdx = zeros(size(allPeaks,1),1);
% find spikes
count = 1;
while count <=size(allPeaks,1)
    k = count+1;
    while (k<=size(allPeaks,1)) && ( (allPeaks(k,1) - allPeaks(count,1)) <= 25)
        k=k+1;
    end
    [mx mxidx] = max(allPeaks(count:k-1,2));
    keepIdx(count+mxidx-1) = 1;
    count=k;
end
allPeaks = allPeaks(keepIdx==1,:);

%% GDF single spikes
gdf_single = allPeaks(:, [3 1]);
R_single = mysort.spiketrain.alignGDFs(tgdf, gdf_single, 5, 3, 25);
mysort.plot.printEvaluationTable(R_single)
% [tgdf(1:10,:) gdf_single(1:10,:)]

mysort.plot.alignment(R_single, 'restrict2Time', [0 1000])

%% PLOT FPs SINGLE !
def = mysort.util.defs();
fpidx = find(R_single.spikeLabel2{1} == def.fp);
fp_times = R_single.St2{1}(fpidx);
for k=1:min(5,length(fp_times))
    fpt = fp_times(k);
    xrange = fpt-400:fpt+400;
    apidx      = ((fpt-400) < allPeaks(:,1)) & (allPeaks(:,1) < (fpt+400));
    tgdfidx    =  (fpt-400) < tgdf(:,2)   &    tgdf(:,2) < (fpt+400);
    gdf_singleidx =  (fpt-400) < gdf_single(:,2)   &    gdf_single(:,2) < (fpt+400);
    
    figure
    plot(xrange,Y(:,xrange)')
    hold on
%     for i=1:3
%         idx = ((fpt-400) < peaksOvp{i}(:,1)) & (peaksOvp{i}(:,1) < (fpt+400));
%         plot(peaksOvp{i}(idx,1), peaksOvp{i}(idx,2), 'o', 'markersize', 14, 'color', cols{i});
%     %     mysort.plot.vectorColor(peaksOvp{i}(:,3)))
%     end

    plot(allPeaks(apidx,1), allPeaks(apidx,2), 'cx', 'markersize', 12, 'linewidth', 2)
    
    mysort.plot.spiketrain(tgdf(tgdfidx,:), 'yoffset', -100, 'colorList', cols)
    mysort.plot.spiketrain(gdf_single(gdf_singleidx,:), 'colorList', cols)
    plot(allPeaks(apidx,1), allPeaks(apidx,2), 'cx', 'markersize', 12, 'linewidth', 2)
    mysort.plot.dataCursorMode();
    axis tight
end



