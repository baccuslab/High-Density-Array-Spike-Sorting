

% %% Buffer short
% X = X(:, 1:400000);
% save('C:\LocalData\Michele\marching_square_buffer_short.mat', 'X','srate','channel_x', 'channel_y');


%% Alternative ALLE
% H = mysort.util.hdf5recursiveLoad('C:\LocalData\Alle\AnalyseAlleCompleteHDF5\Preprocessing\buf0003.h5');
% X1 = H.RESi.r(1,1).X;
% X2 = H.RESi.r(1,2).X;
% save('C:\LocalData\Alle\temp_prep_buf0003h5.mat', 'X1', 'X2');

%% Load data
% % 
% clear all
% close all
% load('C:\LocalData\Michele\marching_square_buffer_short.mat', 'X','srate','channel_x', 'channel_y');
% k1 = 1;
% k2 = 4;

% Alternative
% load('C:\LocalData\Alle\temp_prep_buf0001h5.mat', 'X');
% k1 = 1;
% k2 = 2;

% Alternative
% clear all
% close all
% load('C:\LocalData\Alle\temp_prep_buf0003h5.mat', 'X1');
% X = X1; clear X1;
% k1 = 3;
% k2 = 4;

% Alternative
% clear all
% close all
% load('C:\LocalData\Alle\temp_prep_buf0003h5.mat', 'X1', 'X2');
% shift = 1000;
% X = [X1(3,1:end-shift); X2(3,1+shift:end)]; clear X1 X2;
% k1 = 1;
% k2 = 2;

%% Get sigma
Tf = 100;
thr1 = 4.5;
[smad1, spikeepochs1] = ana.douglas.estimateSigma(X(k1,:), Tf, thr1);
noiseepochs1 = mysort.epoch.flip(spikeepochs1, size(X,2));
[smad2, spikeepochs2] = ana.douglas.estimateSigma(X(k2,:), Tf, thr1);
noiseepochs2 = mysort.epoch.flip(spikeepochs2, size(X,2));

spikeepochs1_idx = mysort.epoch.toIdx(spikeepochs1);
spikeepochs2_idx = mysort.epoch.toIdx(spikeepochs2);
noiseepochs1_idx = mysort.epoch.toIdx(noiseepochs1);
noiseepochs2_idx = mysort.epoch.toIdx(noiseepochs2);
commonnoiseepochs = mysort.epoch.intersect(noiseepochs1,noiseepochs2);
commonnoiseepochs_idx = mysort.epoch.toIdx(commonnoiseepochs);
unionspikeepochs  = mysort.epoch.merge( [spikeepochs1; spikeepochs2]);
unionspikeepochs_idx = mysort.epoch.toIdx(unionspikeepochs);
%% show test of data
fig = mysort.plot.figure();
h1 = subplot(2,1,1);
plot(X(k1,:), 'k'); hold on

thr2 = thr1;
plot([0 size(X,2)], thr2*[smad1 smad1], ':g', 'linewidth', 2)
plot([0 size(X,2)],-thr2*[smad1 smad1], ':g', 'linewidth', 2)
mysort.plot.epochs(unionspikeepochs, 0, 'm', 'linewidth', 6);
mysort.plot.epochs(spikeepochs1, 0, 'c', 'linewidth', 4);
mysort.plot.epochs(noiseepochs1, 0, 'r', 'linewidth', 4);
mysort.plot.epochs(commonnoiseepochs, 0, 'k', 'linewidth', 2);


h2 = subplot(2,1,2);
plot(X(k2,:), 'k'); hold on

plot([0 size(X,2)], thr2*[smad2 smad2], ':g', 'linewidth', 2)
plot([0 size(X,2)],-thr2*[smad2 smad2], ':g', 'linewidth', 2)
mysort.plot.epochs(unionspikeepochs, 0, 'm', 'linewidth', 6);
mysort.plot.epochs(spikeepochs2, 0, 'c', 'linewidth', 4);
mysort.plot.epochs(noiseepochs2, 0, 'r', 'linewidth', 4);
mysort.plot.epochs(commonnoiseepochs, 0, 'g', 'linewidth', 2);

linkaxes([h1 h2],'xy');
%% Estimate autocov
maxLag = 50;
L = size(X,2);
r1 = xcorr(X(k1,:), maxLag, 'biased');
r2 = xcorr(X(k2,:), maxLag, 'biased');
xd = xcorr(X(k1,:), X(k2,:), maxLag, 'biased');

[c1 Ln1] = mysort.util.xcorr_in_epochs(X(k1,:), noiseepochs1, maxLag, maxLag);
[c2 Ln2] = mysort.util.xcorr_in_epochs(X(k2,:), noiseepochs2, maxLag, maxLag);
[xc Lxn] = mysort.util.xcorr_in_epochs(X([k1 k2],:),commonnoiseepochs, maxLag, maxLag);

[sp1 Lsp1] = mysort.util.xcorr_in_epochs(X(k1,:),spikeepochs1, maxLag, maxLag);
[sp2 Lsp2] = mysort.util.xcorr_in_epochs(X(k2,:),spikeepochs2, maxLag, maxLag);
[xsp Lxsp] = mysort.util.xcorr_in_epochs(X([k1 k2],:),unionspikeepochs, maxLag, maxLag);


%% Plot covs
fig = mysort.plot.figure();
xr = -maxLag:maxLag;
plot(xr, r1, 'k:', 'linewidth', 2); hold on
plot(xr, r2, 'c:', 'linewidth', 2);
plot(xr, xc, 'm:', 'linewidth', 2);

plot(xr, c1, 'k', 'linewidth', 2); 
plot(xr, c2, 'c', 'linewidth', 2);
plot(xr, xd, 'm', 'linewidth', 2);

plot(xr, sp1, 'k--', 'linewidth', 2); 
plot(xr, sp2, 'c--', 'linewidth', 2);
plot(xr, xsp, 'm--', 'linewidth', 2);
legend('datacov1', 'datacov2', 'data_xcov', 'noisecov1', 'noisecov2', 'noise_xcov', 'spcov1', 'spcov2', 'sp_xcov');

%% Test relationships
fig = mysort.plot.figure();
plot(xr, c1*Ln1, 'b');
hold on
plot(xr, sp1*Lsp1, 'r');
plot(xr, r1*L, 'g');
plot(xr, c1*Ln1 + sp1*Lsp1, 'k');
residual = c1*Ln1 + sp1*Lsp1 - r1*L;
plot(xr, residual, 'k:');

legend('N', 'SP', 'R', 'N+SP', 'err');


%%
tau = 1;
x = [X(k1,:)/smad1; X(k2,:)/smad2];
lims = 15*[-1 1];
%% Plot instantanen template space between noise in both channels
fig = mysort.plot.figure('width', 740, 'height', 800);
h = mysort.plot.subplot2([2,2]);

[Ck1k1_1 Hk1k1_1] = mysort.noise.plot2(h(1), lims, x, noiseepochs1_idx, spikeepochs1_idx, tau, 1, 1, true);
[Ck2k2_1 Hk2k2_1] = mysort.noise.plot2(h(3), lims, x, noiseepochs2_idx, spikeepochs2_idx, tau, 2, 2, false);
[Ck1k2_0 Hk1k2_0] = mysort.noise.plot2(h(2), lims, x, commonnoiseepochs_idx, unionspikeepochs_idx,   0, 1, 2, false);
[Ck1k2_1 Hk1k2_1] = mysort.noise.plot2(h(4), lims, x, commonnoiseepochs_idx, unionspikeepochs_idx, tau, 1, 2, false);
 


%% OLD PLOT
% fig = mysort.plot.figure('width', 1000, 'height', 600);
% h1 = subplot(2,4,1);
% Ck1k1_1 = mysort.noise.plot(h1, lims, x, noiseepochs1_idx, tau, 1, 1, 'k.', 'noise');
% h2 = subplot(2,4,2);
% Hk1k1_1 = mysort.noise.plot(h2, lims, x, spikeepochs1_idx, tau, 1, 1, 'r.', 'spikes');
% linkaxes([h1 h2], 'xy');
% 
% h1 = subplot(2,4,5);
% Ck2k2_1 = mysort.noise.plot(h1, lims, x, noiseepochs2_idx, tau, 2, 2);
% h2 = subplot(2,4,6);
% Hk2k2_1 = mysort.noise.plot(h2, lims, x, spikeepochs2_idx, tau, 2, 2, 'r.');
% linkaxes([h1 h2], 'xy');
% 
% h1 = subplot(2,4,3);
% Ck1k2_0 = mysort.noise.plot(h1, lims, x, commonnoiseepochs_idx, 0, 1, 2);
% h2 = subplot(2,4,4);
% Hk1k2_0 = mysort.noise.plot(h2, lims, x, unionspikeepochs_idx, 0, 1, 2, 'r.');
% linkaxes([h1 h2], 'xy');
% 
% h1 = subplot(2,4,7);
% Ck1k2_1 = mysort.noise.plot(h1, lims, x, commonnoiseepochs_idx, tau, 1, 2);
% h2 = subplot(2,4,8);
% Hk1k2_1 = mysort.noise.plot(h2, lims, x, unionspikeepochs_idx, tau, 1, 2, 'r.');
% linkaxes([h1 h2], 'xy');

