% Make groupings of electrodes to be sorted independently

electrodePositions = 500*rand(1000,2);

tic
[groupsidx nGroupsPerElectrode] = mysort.mea.constructLocalElectrodeGroups(electrodePositions(:,1), electrodePositions(:,2));
toc