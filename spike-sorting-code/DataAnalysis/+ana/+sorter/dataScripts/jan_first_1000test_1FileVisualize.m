pd = pdefs();
expPath = fullfile(pd.serverData, '..', 'Mea1k', 'shared', '140821', 'data', 'hamster_DS_4');

E = struct();
E.expName = 'hamster_DS_4';
E.mapFile   = fullfile(pd.serverData, '..', 'Mea1k', 'shared', '140821', 'map.mat');
E.dataFiles = {fullfile(expPath, '0000.raw.h5')};
E.sortingOutTempFolder = fullfile(pd.localDataOnNetwork, 'Mea1k', 'shared', '140821TmpOut_DS_4_All_ResultsV3');
E.sortingName = 'Sorting1';
E.finalResultOutFolder = fullfile(pd.localDataOnNetwork, 'Mea1k', 'shared', '140821_hamster_DS_4_ResultsV3');

HDS = mysort.mea.HDSorter(E);
HDS.prefilterFiles()
FL = fullfile(E.sortingOutTempFolder, HDS.tempPrefilteredFiles{1});
DS = mysort.mea.CMOSMEA(FL, 'useFilter', 0, 'name', 'PREFILT'); 

resultFile = 'S:\group\hierlemann\Temp\FelixFranke\LocalData\Mea1k\shared\140821_hamster_DS_4_Results\Sorting1_results.mat';
R = load(resultFile);

M = load(E.mapFile);

%%


% %%
% SDA = mysort.plot.SliderDataAxes({DS});
% SDA.setXLim([0 50]);
% 
% %%
% mysort.plot.templates2D(R.T_merged/5, ep, 4)

%%
ep = DS.MultiElectrode.electrodePositions;
% figure, hist(ep(:,2), 200)
% blocSepY = 1200;
% blocSepX = 2200; 
% dimIdx = 1;
% dX = min(ep(ep(:,dimIdx)>blocSepX,1)) - min(ep(ep(:,dimIdx)<blocSepX,1));
% dimIdx = 2;
% dY = min(ep(ep(:,dimIdx)>blocSepY,1)) - min(ep(ep(:,dimIdx)<blocSepY,1));
% 
% rightElBlockIdx = find(ep(:,dimIdx)>blocSepX);
% 
% ep(rightElBlockIdx, :) = ep(rightElBlockIdx, :) + repmat([-dX 0], length(rightElBlockIdx),1 );

ux = unique(ep(:,1));
uy = unique(ep(:,2));

[~, epix] = ismember(ep(:,1), ux);
[~, epiy] = ismember(ep(:,2), uy);



% Plot the gridded data as a mesh and the scattered data as dots.

KERN =  [-1 -1 -1 -1 -1 -1 -1
         -1  0  0  0  0  0 -1
         -1  0  2  2  2  0 -1
         -1  0  2  3  2  0 -1 
         -1  0  2  2  2  0 -1
         -1  0  0  0  0  0 -1
         -1 -1 -1 -1 -1 -1 -1];
KERN = KERN/norm(KERN(:));
 
siz = [length(ux), length(uy)];
idx  = sub2ind(siz, epix, epiy);
subidx = ind2sub(siz, 1:(length(ux)*length(uy)));
nidx = sub2ind(siz, setdiff(subidx, idx));
U = unique(R.gdf_merged(:,1));

range = [-1 -1];
[mi, ma] = mysort.wf.tMinMaxPerTemplate(R.T_merged);
P2P = ma-mi;
%%
fh = mysort.plot.figure([1200 800]);
ah = subplot(4,1,1);

ah(2) = subplot(4,1,2:4);
axis equal
colormap(gray)
for i=5000:1:size(DS,1)
    x = DS(i,:)/10^4;
    
    X = ones(length(ux), length(uy))*mean(x);
    X(idx) = x;
    
    if i<range(1) || i > range(2)
        cla(ah(1))
        range = [i i+999];
        plot(ah(1), range(1):range(2), DS(i:i+999,:));
        set(ah(1), 'nextplot', 'add');
        
    end
        
    
%     X = conv2(X, KERN, 'same');
    X(idx) = X(idx)-mean(x);
%     X(nidx) = 0;
    
    cla(ah(2)); 
    try
        delete(lh)
    end
    try
        for k=1:length(P.LH)
            delete(P.LH{k})
        end
    end
%     [out, cb] = mysort.plot.imagesc(X, [-300 300]); 
% (1:22,1:22)
    out = mysort.plot.imagesc(ah(2), X', [-1500 1500]);
    set(out, 'xdata', [min(ux) max(ux)], 'ydata', [min(uy) max(uy)]);
%     closeSpikes = find(abs(R.gdf_merged(:,2)-i) < 8);
%     if ~isempty(closeSpikes)
%         [~, closeIDs] = ismember(R.gdf_merged(closeSpikes,1), U);
%         P = mysort.plot.templates2D(R.T_merged(:,:,closeIDs)/5, ep, 20, 50, 'plotLegend', 0, 'ah', ah(2));
%     end
    axis equal
    lh = mysort.plot.verticalLines(ah(1), i);
%     SDA.setSliderPosition(i, 1)
%     SDA.bAvoidReentrant = false;
%     set(cb, 'YLim' , [-300 300]);
%     set(gca, 'zlim', [-300 300]);
%     set(gca, 'xlim', [min(ep(:,1)) max(ep(:,1))]);
%     set(gca, 'ylim', [min(ep(:,2)) max(ep(:,2))]);
    drawnow
    pause(.01)
end

%%
figure
axes
for i= 1:size(R.T_merged,3)
    if ismember(i, [2 4 5 7 9 10 ])
        continue
    end
    for kk=1:1
        for k=1:size(R.T_merged,1)
            x = squeeze(R.T_merged(k, :, i));
    %         X = ones(length(ux), length(uy))*mean(x);
            X = zeros(length(ux), length(uy))*mean(x);
            X(idx) = x;     
    %         X(idx) = X(idx)-mean(x);
            imagesc(X, [-400 400]);
            title(sprintf('Neuron ID: %d  Framenumber: %d', i, k));
            drawnow
            pause(.06)
        end
    end
end
%%
count = 0;
ah = [];
for i= 1:size(R.T_merged,3)
    if ismember(i, [2 4 5 7 9 10 ])
        continue
    end
    if count == length(ah)
        mysort.plot.figure([1500 1000])
        ah = mysort.plot.subplots([3 5], 'spacerX', .01, 'spacerY', .01);
        count = 0;
    end
    count = count+1;
    x = zeros(length(ux), length(uy));
    x(idx) = P2P(i,:);
    imagesc(x, 'parent', ah(count));
    set(ah(count), 'xticklabel', [], 'yticklabel', [])
    axis(ah(count), 'ij') %,  'equal'
end