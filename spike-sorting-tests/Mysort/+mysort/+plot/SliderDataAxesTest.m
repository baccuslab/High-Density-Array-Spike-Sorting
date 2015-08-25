%%
X1 = randn(10000,10);
mysort.plot.SliderDataAxes(X1);

%%
X1 = randn(10000,10);
X2 = randn(10000,10);
mysort.plot.SliderDataAxes({X1, X2});

%%
X1 = randn(10000,10);
samplesPerSecond1 = 32000;
ds1 = mysort.ds.Matrix(X1, samplesPerSecond1);
mysort.plot.SliderDataAxes(ds1, 'timeIn', 'ms');

%%
X1 = randn(10000,10);
samplesPerSecond1 = 32000;
ds1 = mysort.ds.Matrix(X1, samplesPerSecond1);
X2 = randn(10000,10);
samplesPerSecond2 = 1000;
ds2 = mysort.ds.Matrix(X2, samplesPerSecond2);
mysort.plot.SliderDataAxes({ds2 ds1}, 'timeIn', 'ms');


%%
close all
clear all

X1 = randn(32000,10);
samplesPerSecond1 = 32000;
ds1 = mysort.ds.Matrix(X1, samplesPerSecond1);
T = [];
for i=1:10
    T(:,i,1) = i*sin(0:.2:3);
    T(:,i,2) = i*cos(0:.2:3);
end
cutLeft = 3;
gdf = [ones(10,1)*1 linspace(100, 9000, 10)';
       ones(10,1)*2 linspace(100, 9000, 10)'];
Sorting = mysort.spiketrain.SpikeSortingContainer('Sorting1', gdf, 'templateWfs', T, 'templateCutLeft', cutLeft);
ds1.addSpikeSorting(Sorting);

X2 = randn(3000,10);
samplesPerSecond2 = 1000;
ds2 = mysort.ds.Matrix(X2, samplesPerSecond2);
mysort.plot.SliderDataAxes({ds2 ds1}, 'timeIn', 'ms', 'dataColors', {[], 'k'}, 'plotSortings', 1);