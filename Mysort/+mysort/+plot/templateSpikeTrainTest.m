T = [];
T(:,:,1) = [1:10; 10:-1:1]';
T(:,:,2) = [10:-1:1; 1:10]';

gdf = [3 15
       9 33
       3 50
       3 55
       9 57];

tM = 1;
chSp = 10;
plot_mode = 'normal';
cutLeft = 2;
unitNames = [3 9];

figure;
mysort.plot.templateSpikeTrain(T, gdf,...
        'timeMultiplicator', tM, ...
        'channelSpacer', chSp, ...
        'mode', plot_mode,...
        'cutLeft', cutLeft,...
        'T_gdf_idx2id', unitNames);