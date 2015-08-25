%% INIT TEST FILE
cd('/net/bs-filesvr01/export/group/hierlemann/Temp/FelixFranke/LocalData/Michele')
ffile = 'gratings.stream.ntk';

outFile1 = 'test_gratings11Filtered.stream.h5';

%% FIRST CONVERT DATA TO H5
if ~exist(outFile1, 'file');
    mysort.mea.convertNTK2HDF(ffile, 'outFile', outFile1, 'prefilter', 1);
end

%% INIT DATASOURSE
ds_mea = mysort.mea.CMOSMEA(outFile1, 'useFilter', 0, 'name', 'DATA');

%% START SORTING
outPath = 'TestSorting';
[gdf_merged T_merged localSorting localSortingID] = mysort.groupedSorting.startHDSorting(ds_mea, outPath);

%% PLOT ALL TEMPLATES ON TOP OF EACH OTHER
nShowElectrodes = 5;
showOnlyElectrodesWithMinAbsAmplitudeOf = [];
mysort.plot.templates2D(T_merged, ds_mea.MultiElectrode.electrodePositions,...
    nShowElectrodes,...
    showOnlyElectrodesWithMinAbsAmplitudeOf);

%% PLOT EACH TEMPLATE INDIVIDUALLY 
nShowElectrodes = 25;
mysort.plot.templates2D(T_merged, ds_mea.MultiElectrode.electrodePositions,...
    nShowElectrodes,...
    showOnlyElectrodesWithMinAbsAmplitudeOf,...
    'stacked', 0);