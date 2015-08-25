%% INIT FELIX
% matlabpool(12);
pd = pdefs();
% make sure you have the path of a folder with all .ntk files WITH A SINGLE
% COMMON configuration in the variable "dpath"
expNames = {'20140826'};

        

for kk=1:length(expNames)
expName = expNames{kk};
sortingName = 'Sorting1';
depath = fullfile(pd.networkTempShare, 'Tom', 'Sorting');
outPath = fullfile(depath, [expName 'Out']);
if ~exist(outPath, 'file'); mkdir(outPath); end
confListFile = fullfile(depath, [expName '.mat']);
c = load(confListFile);
configNTKLists = c.configNTKLists;

%% Now we CONVERT all ntk files to H5. Doing this we also prefilter the data
h5FileLists = {};
for i=1:length(configNTKLists)
    for k=1:length(configNTKLists{i})
        fn = configNTKLists{i}{k};
        [ffolder ffname fsuffix] = fileparts(fn);
        fnpf = fullfile(outPath, [ffname fsuffix '_prefilter.h5']);
        h5FileLists{i}{k} = fnpf;
        bConvertFile = true;
        if exist(fnpf, 'file')
            try
                bIsInProcess = hdf5read(fnpf, '/bFileIsInProcess');
%                 bIsInProcess = 1;
                if bIsInProcess
                    % This means the file is still in process. We assume no
                    % other user/matlab instance is writing so it is
                    % probably a left over from some other run and we
                    % delete it.
                    delete(fnpf); 
                else
                    % The file was already properly processed, we dont need
                    % to do anything.
                    bConvertFile = false;
                end
            catch
                delete(fnpf);
            end
        end
        if bConvertFile
            % Only convert file if it was non existent or deleted
            mysort.mea.convertNTK2HDF(fn, 'prefilter', 1, 'outFile', fnpf);
        end
    end
end

%%
if 0
for i=1:length(h5FileLists)
    sortOutPath = fullfile(outPath, ['Config' num2str(i)], sortingName);
    ana.michele.spikeSortingStartScripts.run(h5FileLists{i}, sortOutPath, sortingName); 
end
end

%% EXPORT TO SINGLE FOLDER
exportPath = fullfile(outPath, 'Results');
if ~exist(exportPath, 'file')
    mkdir(exportPath)
end
R = {};
for i=1:length(h5FileLists)
    sortOutPath = fullfile(outPath, ['Config' num2str(i)], sortingName);
    sourceFile = fullfile(sortOutPath, [sortingName '_results.mat'])
    assert(exist(sourceFile, 'file')>0, 'File not found!');
    D = load(sourceFile);
    
    % get individual files' length in samples
    compoundMea = mysort.mea.compoundMea(h5FileLists{i}, 'useFilter', 0, 'name', 'PREFILT');
    L = compoundMea.X.getAllSessionsLength();
    
    % break gdf apart
    start_end_times = [0 cumsum(L)];
    mgdf = mysort.spiketrain.gdf2multiSessionGdf(D.gdf_merged, start_end_times);
    for k=1:length(L)
        gdf = mgdf(mgdf(:,3)==k,[1:2]);
        R{k, i} = gdf;
    end
end
save(fullfile(depath, [expName '_resultsForMichele.mat']), 'R');
end
% 
% 

%% 
outPath = fullfile(depath, [expName 'Out']);

configFolderPathList = {};
configNTKLists2 = {};
for i=1:length(h5FileLists)
    for j=1:length(h5FileLists{i})
        fn = configNTKLists{i}{j};
        [ffolder ffname fsuffix] = fileparts(fn);
        fnpf = [ffname fsuffix];
        configNTKLists2{i}{j} = fullfile(outPath, fnpf);
    end
    configFolderPathList{i} = fullfile(outPath, 'Config1');
end

%%
configIdx = 1;

[SSC S compoundMea G] = ana.toni.inspectSorting(sortingName, configFolderPathList{configIdx}, configNTKLists2{configIdx});
unitsFound = SSC.unitNames;
cutLeft = 10;
cutLength = 45; 

%%
showUnitsIdx = [2 12 20];
[~, showUnitsIdx] = ismember([902], unitsFound); % 903 904 905 906 907], unitsFound);
for i=1:length(showUnitsIdx)
    nElectrodes = size(compoundMea, 2);
    uIdx = showUnitsIdx(i);
    uID = unitsFound(uIdx);
    uElGroup = (uID - mod(uID,100))/100;
    spikes = SSC.getWaveforms4UnitIdx(uIdx, cutLeft, cutLength);
    tspikes = mysort.wf.v2t(spikes, nElectrodes);
    EP = compoundMea.MultiElectrode.electrodePositions;
    groupElIdx = G.groupsidx{uElGroup};
    figure
    mysort.plot.waveforms2D(tspikes, EP, 'plotMedian', 1, 'plotElNumbers', compoundMea.MultiElectrode.electrodeNumbers);
    hold on
    plot(EP(groupElIdx,1), EP(groupElIdx,2), 'rx', 'markersize', 14, 'linewidth', 3)
    mysort.plot.figureTitle(['UnitIdx: ' num2str(uIdx) ' UnitName:' num2str(unitsFound(uIdx))]);
end

%%
selectedElectrodeNumbers = [240 342 241 137 139 242 344];
[~, selectedElectrodeIndices] = ismember(selectedElectrodeNumbers, compoundMea.MultiElectrode.electrodeNumbers);
compoundMea.restrictToChannels(selectedElectrodeIndices);
% compoundMea.restrictToChannels([]);
%%
mysort.plot.SliderDataAxes(compoundMea, 'plotSortings', 0, 'channelSpacer', 50, 'timeIn', 's');


%%
X = compoundMea(:, :);
Y = sum(X,2);
%%
figure, plot(Y)
%%
tic
[pks, TS] = findpeaks(-Y,'MinPeakHeight', 175, 'MinPeakDistance', 50);
toc
%%
compoundMea.restrictToChannels([]);
spikes = compoundMea.getWaveform(TS, 25, 75);
tSpikes = mysort.wf.v2t(spikes, size(compoundMea,2));

%%
S1 = squeeze(tSpikes(:, selectedElectrodeIndices(1), :));
figure
plot(S1)

%%
A1a = min(S1,[],1)';
A1b = max(S1,[],1)';
removeIdx = ((A1b < 50 | A1a > -50) & TS < 3*10^6)  |  ((A1b < 50 | A1a > -50) & TS > 6.1*10^6);
TS_ = TS(~removeIdx);
figure
plot(TS_, [A1a(~removeIdx) A1b(~removeIdx)], '.')

%%
S1_ = S1(:,~removeIdx);
%%
figure
surf(S1_(:,1:10:end));
%%
L = 200;
S1f = zeros(size(S1_));
for i=1:size(S1_,1)
    S1f(i,:) = conv(S1_(i,:), ones(L,1)/L, 'same');
end

%%
figure
surf(S1f(:,1:50:end));

%%
temperatureSteps = 400000;
tSpikes_ = tSpikes(:,:, ~removeIdx);
trange = temperatureSteps:temperatureSteps:size(compoundMea,1);
T = [];
for i=1:length(trange)
    idx = TS_ > trange(i)-temperatureSteps & TS_ < trange(i);
    if sum(idx) > 1
        T(:,:,i) = mean(tSpikes_(:,:,idx), 3);
    end
end

%%
temperatures = [20:-.5:15 15.5:.5:20];
miT = min(temperatures);
maT = max(temperatures);
rt = maT-miT;
figure
ah = axes;
hold on
for i=1:length(temperatures);
    t = temperatures(i);
    ti = (t-miT)/rt;
    col = [1 0 0]*ti + (1-ti)*[0 0 1];
    mysort.plot.waveforms2D(T(:,:,i), EP, 'AxesHandle', ah, 'plotArgs', {'linewidth', 2, 'color', col}) %, 'plotMedian', 1, 'plotElNumbers', compoundMea.MultiElectrode.electrodeNumbers);
end
    
%%
ts = D.gdf_merged(D.gdf_merged(:,1)==902,2);
figure
plot(ts(2:end), diff(ts), '.')

%%
figure
plot(ts/20000, ones(1, length(ts)), '.')


