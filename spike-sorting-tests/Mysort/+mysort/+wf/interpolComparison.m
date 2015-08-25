% This example shows that the effect of sampling jitter can be best
% cancelled using the resample logic of matlab. neither interpolating with
% splines nor interpolating with sinc filters gives in this case a good
% outcome. The high frequency parts, that are over the downsampled nyquest 
% frequency in the simulated test are of course lost for all methods.

%% create test data
ds = 100;
maxShift = ds;
y = 3*[zeros(1, maxShift) 0:.01:10 10.5 10-.01:-.01:0 zeros(1,maxShift)];
xorg = 1:length(y);
T = y(1:end-maxShift+1);
Y = resample(T, 1, ds);
Y(maxShift, end) = 0;
T(maxShift, end) = 0;
for i=2:maxShift
    T(i,:) = y(i:end-maxShift+i);
    Y(i,:) = resample(T(i,:), 1, ds);
end
xds = resample(xorg, 1, ds);
shifts = (1:maxShift)/ds;
shifts = -(shifts'-median(shifts));

%% re-upsample testdata
Z = resample(Y(1,:), ds, 1);
Z(maxShift, end) = 0;
for i=2:maxShift
    Z(i,:) = resample(Y(i,:), ds, 1);
end

%% interpolate data
Q = [];
for i=1:maxShift
    ipf = mysort.util.interpolfun(Y(i,:));
    Q(i,:) = ipf(1:.01:size(Y,2));
end

%% sinc interpolate data
R = [];
for i=1:maxShift
    sf = mysort.wf.mSincfun(Y(i,:));
    R(i,:) = sf(1:.01:size(Y,2));
end

%% PLOT
figure;p=0;nP=5;ah=zeros(nP,1);

p=p+1; ah(p) = subplot(nP,1,p);
plot(T', '.-');
title('raw data');

p=p+1; ah(p) = subplot(nP,1,p);
plot(Y', '.-');
title('downsampled data');

p=p+1; ah(p) = subplot(nP,1,p);
plot(Z', '.-');
title('re-upsampled data');

p=p+1; ah(p) = subplot(nP,1,p);
plot(Q', '.-');
title('re-interpolated data');

p=p+1; ah(p) = subplot(nP,1,p);
plot(R', '.-');
title('re-sincinterpolated data');

linkaxes(ah, 'y');
set(ah, 'ylim', [29.5 30.1]);

%%
% min_x_i = zeros(size(Y,1),1);
% min_y_i = min_x_i;
% min_x = min_x_i;
% min_y = min_x_i;
% for i=1:size(Y,1)
%     [min_x_i(i) min_y_i(i) min_x(i) min_y(i)] = mysort.util.findMinInterp(-Y(i,:));
% end
% min_x_i = min_x_i-median(min_x_i);
% % [tau2 Ya2] = mysort.util.alignWaveformsInterpolatedMean(Y, 1);
% % subplot(3,1,3)
% % plot(Ya2');
% % title('aligned on max');
% fprintf('Align on Interp Max error: %f\n', sum(abs(min_x_i-shifts)));
% 
% %%
% figure;
% plot(shifts)
% hold on
% plot(min_x_i, 'g');
% plot(-tau1, 'r');
% legend('simulated shifts','shifts estimated by findMinInterp','shifts estimated by vAlignOnMax');
% 
% 
% % findMinInterp