%% INIT FELIX
pd = pdefs();
% make sure you have the path of a folder with all .ntk files WITH A SINGLE
% COMMON configuration in the variable "confPath"

% x = EL.map.mposx * 17.5;
% y = EL.map.mposy * 17.5;
% 
% figure;scatter( EL.map.mposx , EL.map.mposy )
% axis ij
% axis equal

% matlabpool(4)

expPath = fullfile(pd.serverData, '..', 'Mea1k', 'shared', '140821', 'data', 'hamster_DS_4');

flist = dir(fullfile(expPath, '*.h5'));
FL = {};
for i=1:length(flist)
    FL{i} = fullfile(expPath, flist(i).name);
end

P = struct();
P.expName = 'hamster_DS_4All';
P.mapFile   = fullfile(pd.serverData, '..', 'Mea1k', 'shared', '140821', 'map.mat');
P.dataFiles = FL;
P.sortingOutTempFolder = fullfile(pd.localData, 'Mea1k', 'shared', '140821TmpOut_DS_4_All_Results2');
P.sortingName = 'Sorting1';
P.finalResultOutFolder = fullfile(pd.localData, 'Mea1k', 'shared', '140821_hamster_DS_4_All_Results2');

M = load(P.mapFile);

HDS = mysort.mea.HDSorter(P);
% matlabpool close
% matlabpool(12)
disp('Prefiltering')
HDS.prefilterFiles()
disp('Done.')
% oneFile = fullfile(P.sortingOutTempFolder, HDS.tempPrefilteredFiles{1});
% DS = mysort.mea.CMOSMEA(oneFile, 'useFilter', 0, 'name', 'PREFILT'); 
% ep = DS.MultiElectrode.electrodePositions;

%%
% matlabpool close
% matlabpool(8)
t = tic;
HDS.runSorting();
toc(t)

return

%%
% map = M.map;
% isConnected = map.elNo>0;
% nC_effective = length(find(isConnected));
% for i=3:3 %4:length(HDS.tempPrefilteredFullFiles)
%     i   
%     P.outFile = HDS.tempPrefilteredFullFiles{i};
%     idx = isConnected==1;
%     idx = idx(:);
%     x = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/channel_nr', [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
%     x(1,1:nC_effective) = map.elNo(idx)';
%     clear x
%     x = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/channel_posx', [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
%     x(1,1:nC_effective) = map.mposx(idx)';
%     clear x
%     x = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/channel_posy', [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
%     x(1,1:nC_effective) = map.mposy(idx)';
%     clear x
%     x = mysort.h5.createVariableAndOrFile(P.outFile, '/Sessions/Session0/channel_connected', [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
%     x(1,1:nC_effective) = ones(1,nC_effective);
%     clear x
% end
% return 
%%
% DSFull = mysort.mea.compoundMea(HDS.tempPrefilteredFullFiles, 'useFilter', 0, 'name', 'PREFILT'); 

%%

% %%
% figure
% plot(ep(:,1), ep(:,2), 'x')
% 
% %%
% 
% 
% %%
% figure
% for i=1:size(DS,1)
%     x = DS(i,:);
%     clf 
%     surf(ep(:,1), ep(:,2), x);
%     pause(.2)
% end
% %%
% figure
% plot(DS(1:2000,:))

%%
ux = unique(ep(:,1));
uy = unique(ep(:,2));

[~, epix] = ismember(ep(:,1), ux);
[~, epiy] = ismember(ep(:,2), uy);

% Plot the gridded data as a mesh and the scattered data as dots.

k = [1 1 1
     1 1 1 
     1 1 1]/9;
 
siz = [length(ux), length(uy)];
idx  = sub2ind(siz, epix, epiy);
subidx = ind2sub(siz, 1:(length(ux)*length(uy)));
nidx = sub2ind(siz, setdiff(subidx, idx));

figure
colormap(gray)
for i=1:size(DS,1)
    x = DS(i,:);
    
    X = ones(length(ux), length(uy))*mean(DS(i,:));
    X(idx) = DS(i,:);
    X = conv2(X, k, 'same');
    X(idx) = X(idx)-mean(DS(i,:));
    X(nidx) = 0;
    
    clf 
%     [out, cb] = mysort.plot.imagesc(X, [-300 300]);
    mysort.plot.imagesc(X, [-1500 1500]);
%     set(cb, 'YLim' , [-300 300]);
%     set(gca, 'zlim', [-300 300]);
%     set(gca, 'xlim', [min(ep(:,1)) max(ep(:,1))]);
%     set(gca, 'ylim', [min(ep(:,2)) max(ep(:,2))]);
    pause(.01)
end
%%
figure
plot(DS(1:20000,:))

%%
mysort.plot.SliderDataAxes({DS})



% M = 