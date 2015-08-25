
clear all
close all
%%
load('C:\LocalData\Michele\marching_square_buffer_short.mat', 'X','srate','channel_x', 'channel_y');
maxDist = 300;
% [xcovs D autocovs autocov_norms NE smad] = mysort.mea.noise_covariance_over_distance(X,channel_x,channel_y, 'maxDist', maxDist);

michele = load('buffer_michele', 'xcovs', 'D', 'autocovs', 'autocov_norms', 'NE', 'smad');


%%
[Y srate x y] = ana.douglas.getPreprocessed('BlockEasy');
fig = mysort.plot.figure;
mysort.mea.plotArray(x,y);
% [xcovs D autocovs autocov_norms NE smad] = mysort.mea.noise_covariance_over_distance(Y,x,y, 'maxDist', maxDist);
% save('buffer_doug', 'xcovs', 'D', 'autocovs', 'autocov_norms', 'NE', 'smad');
douglas = load('buffer_doug', 'xcovs', 'D', 'autocovs', 'autocov_norms', 'NE', 'smad');

%%
mysort.plot.figure('width', 800, 'height', 500);
h(1) = subplot(2,1,1);
set(h(1), 'fontsize',14);
max_xc = max(michele.xcovs, [], 2);
% max_xc = michele.xcovs(:,10);
plot(michele.D, max_xc, '.k');
hold on

max_xc = max(douglas.xcovs, [], 2);
% max_xc = douglas.xcovs(:,10);
plot(douglas.D, max_xc, '.r');
xlabel('distance [mu m]');
ylabel('max(covariance)');
title(sprintf('noise cross channel covariance, maxDist: %d, maxlag: 10',maxDist));
legend('Michele', 'Douglas');
set(h(1), 'xlim', [0 300])

h(2) = subplot(2,1,2);
set(h(2), 'fontsize',14);
plot(-10:10, michele.autocovs', 'k')
hold on
plot(-10:10, douglas.autocovs', 'r')
title('auto covariance functions');
xlabel('timelag');
ylabel('covariance');
mysort.plot.savefig(gcf, 'NoiseCorrDougMichele');