%% simple
y = zeros(1,100);
y(51:55) = [-1 -2 -2.5 -.4 2];
[val, idx] = mysort.util.findMinInterp(y);


%% create test data
ds = 100;
maxShift = ds;
y = 3*[zeros(1, maxShift) 0:.01:10 10-.01:-.01:0 zeros(1,maxShift)];
T = y(1:end-maxShift+1);
Y = resample(T, 1, ds);
Y(maxShift, end) = 0;
T(maxShift, end) = 0;
for i=2:maxShift
    T(i,:) = y(i:end-maxShift+i);
    Y(i,:) = resample(T(i,:), 1, ds);
end
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
    Q(i,:) = ipf(1:.001:size(Y,2));
end

%% PLOT
figure;p=0;nP=4;ah=zeros(nP,1);
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
linkaxes(ah, 'y');
%%
[tau1 Ya1] = mysort.util.alignWaveformsOnMax(Y, 1);
subplot(3,1,2)
plot(Ya1', '.-');
title('aligned on max');
fprintf('Align on Max error: %f\n', sum(abs(tau1+shifts)));

%%
min_x_i = zeros(size(Y,1),1);
min_y_i = min_x_i;
min_x = min_x_i;
min_y = min_x_i;
for i=1:size(Y,1)
    [min_x_i(i) min_y_i(i) min_x(i) min_y(i)] = mysort.util.findMinInterp(-Y(i,:));
end
min_x_i = min_x_i-median(min_x_i);
% [tau2 Ya2] = mysort.util.alignWaveformsInterpolatedMean(Y, 1);
% subplot(3,1,3)
% plot(Ya2');
% title('aligned on max');
fprintf('Align on Interp Max error: %f\n', sum(abs(min_x_i-shifts)));

%%
figure;
plot(shifts)
hold on
plot(min_x_i, 'g');
plot(-tau1, 'r');

% findMinInterp