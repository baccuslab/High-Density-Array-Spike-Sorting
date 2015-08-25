% f = mysort.plot.figure('w', 1200, 'h', 800);
% mysort.plot.savefig(f, 'B', 'png', 0);
% guide('B.fig')
% MYSORT_PATH = '/home/frankef/bel.svn/hima_internal/cmosmea/trunk/matlab/';
addpath(fullfile('..', '..', '..', 'Mysort', 'examples', 'E10_Quiroga2004'));
path_to_simulated_files = E10_quirogaDataPath();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data, use the ground truth information to cut out spikes,
% calculate the noise statistics and the correct templates
% Choose a benchmark and a noise level
benchmarks = {'Easy1', 'Easy2', 'Difficult1', 'Difficult2'};
b = 1; noisestr = '005';
quirogaFile = sprintf('C_%s_noise%s.mat', benchmarks{b}, noisestr);

%print the statistics to that benchmark:
E10_printQuirogaEvaluation(b);

% Set the length of the templates and where they should be cut in respect
% to the spike times stored in the ground truth data
Tf = 71; cutLeft = -10;
GT = E10_preprocessing(path_to_simulated_files, quirogaFile, cutLeft, Tf);
b = 4; noisestr = '02';
quirogaFile = sprintf('C_%s_noise%s.mat', benchmarks{b}, noisestr);
GT2 = E10_preprocessing(path_to_simulated_files, quirogaFile, cutLeft, Tf);

nT = 3;
nS = 100000; %length(GT.X);
nC = 1;
gdf = GT.gdf;
gdf2 = GT2.gdf;
Tf = size(GT.XI,1);

figure;
subplot(2,2,1)
plot(squeeze(GT.XIunaligned))
subplot(2,2,3)
plot(squeeze(GT.XI))
subplot(2,2,2)
plot(squeeze(GT2.XIunaligned))
subplot(2,2,4)
plot(squeeze(GT2.XI))

for t=1:nT
    gdf(gdf(:,1)==t,2) = gdf(gdf(:,1)==t,2) + GT.tau(t);
    gdf2(gdf2(:,1)==t,2) = gdf2(gdf2(:,1)==t,2) + GT2.tau(t);
end
%%
Yeucl  = zeros(nS, nT);
Yconv  = zeros(nS, nT);
Ymatch = zeros(nS, nT);
Ybotm = zeros(nS, nT);
for t=1:3
    T = squeeze(GT.XI(:,1,t));
    Yconv(:,t)  = conv2(GT.X(1:nS)', flipud(T)/norm(T), 'same');
    CC = (GT.C + GT.C(1,1)*eye(Tf))/2;
    F = (T'/CC)';
%     F = (T'*inv(GT.C))';
    Ymatch(:,t) = conv2(GT.X(1:nS)', flipud(F), 'same');
    Ybotm(:,t)  = Ymatch(:,t) - .5*F'*T + log(.0001);
    Ymatch(:,t) = Ymatch(:,t)/norm(F);
    Tf2 = round(Tf/2);
    for s=Tf2:nS-Tf2
        s1 = s-Tf2+1;
        s2 = s1+Tf-1;
        Yeucl(s,t) = norm(GT.X(s1:s2)'-T)/nC;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort the data using the correct spike templates and noise estimates
% GT.C = GT.C + eye(size(GT.C))*GT.C(1,1)/10;
% NE = mysort.util.NoiseEstimator(GT.C, Tf);
% botm = mysort.sorters.BOTM(NE, Tf, GT.templates, 'upsample', 3, 'spikePrior', .01);

% PLOT
f = open('B.fig');

axs = mysort.plot.getAxesWithTagFromFigure(f, 'axes1');
axs(2) = mysort.plot.getAxesWithTagFromFigure(f, 'axes6');
axs(3) = mysort.plot.getAxesWithTagFromFigure(f, 'axes7');
axs(4) = mysort.plot.getAxesWithTagFromFigure(f, 'axes8');
axs(5) = mysort.plot.getAxesWithTagFromFigure(f, 'axes9');
axs(6) = mysort.plot.getAxesWithTagFromFigure(f, 'axes10');
axs(7) = mysort.plot.getAxesWithTagFromFigure(f, 'axes11');
axs(8) = mysort.plot.getAxesWithTagFromFigure(f, 'axes12');
axs(9) = mysort.plot.getAxesWithTagFromFigure(f, 'axes13');
axs(10) = mysort.plot.getAxesWithTagFromFigure(f, 'axes14');

plotIDs = [4 20 8];
pgdf = gdf;
pgdf(:,1) = plotIDs(gdf(:,1));
pgdf2 = gdf2;
pgdf2(:,1) = plotIDs(gdf2(:,1));

cols = {mysort.plot.vectorColor(plotIDs(1)) mysort.plot.vectorColor(plotIDs(2)) mysort.plot.vectorColor(plotIDs(3))};
lw = 2;
ax = axs(1);
plot(ax, GT.X, 'k', 'linewidth', lw)
set(ax, 'nextplot', 'add')
mysort.plot.templateSpikeTrain(GT.XI, pgdf, 'AxesHandle', ax, 'T_gdf_idx2id', plotIDs);
ana.botmpaper.axesConfig1(ax, nC);


ax = axs(2);
plot(ax, Yeucl(:,1), 'color', cols{1}, 'linewidth', lw)
set(ax, 'nextplot', 'add');
plot(ax, Yeucl(:,2), 'color', cols{2}, 'linewidth', lw)
plot(ax, Yeucl(:,3), 'color', cols{3}, 'linewidth', lw)
plot(ax, [1 nS], [0 0], ':k', 'linewidth', lw)
ana.botmpaper.axesConfig1(ax, nC);

ax = axs(3);
plot(ax, Yconv(:,1), 'color', cols{1}, 'linewidth', lw)
set(ax, 'nextplot', 'add');
plot(ax, Yconv(:,2), 'color', cols{2}, 'linewidth', lw)
plot(ax, Yconv(:,3), 'color', cols{3}, 'linewidth', lw)
plot(ax, [1 nS], [0 0], ':k', 'linewidth', lw)
ana.botmpaper.axesConfig1(ax, nC);

ax = axs(4);
plot(ax, Ymatch(:,1), 'color', cols{1}, 'linewidth', lw)
set(ax, 'nextplot', 'add');
plot(ax, Ymatch(:,2), 'color', cols{2}, 'linewidth', lw)
plot(ax, Ymatch(:,3), 'color', cols{3}, 'linewidth', lw)
plot(ax, [1 nS], [0 0], ':k', 'linewidth', lw)
ana.botmpaper.axesConfig1(ax, nC);

ax = axs(5);
plot(ax, Ybotm(:,1), 'color', cols{1}, 'linewidth', lw)
set(ax, 'nextplot', 'add');
plot(ax, Ybotm(:,2), 'color', cols{2}, 'linewidth', lw)
plot(ax, Ybotm(:,3), 'color', cols{3}, 'linewidth', lw)
plot(ax, [1 nS], [0 0], ':k', 'linewidth', lw)
ana.botmpaper.axesConfig1(ax, nC);


Yeucl  = zeros(nS, nT);
Yconv  = zeros(nS, nT);
Ymatch = zeros(nS, nT);
Ybotm = zeros(nS, nT);
for t=1:3
    T = squeeze(GT2.XI(:,1,t));
    Yconv(:,t)  = conv2(GT2.X(1:nS)', flipud(T)/norm(T), 'same');
    CC = (GT2.C + GT2.C(1,1)*eye(Tf))/2;
    F = (T'/CC)';
%     F = (T'*inv(GT.C))';
    Ymatch(:,t) = conv2(GT2.X(1:nS)', flipud(F), 'same');
    Ybotm(:,t)  = Ymatch(:,t) - .5*F'*T + log(.0001);
    Ymatch(:,t) = Ymatch(:,t)/norm(F);
    Tf2 = round(Tf/2);
    for s=Tf2:nS-Tf2
        s1 = s-Tf2+1;
        s2 = s1+Tf-1;
        Yeucl(s,t) = norm(GT2.X(s1:s2)'-T)/nC;
    end
end

ax = axs(6);
plot(ax, GT2.X, 'k', 'linewidth', lw)
set(ax, 'nextplot', 'add')
mysort.plot.templateSpikeTrain(GT2.XI, pgdf2, 'AxesHandle', ax, 'T_gdf_idx2id', plotIDs);
ana.botmpaper.axesConfig1(ax, nC);


ax = axs(7);
plot(ax, Yeucl(:,1), 'color', cols{1}, 'linewidth', lw)
set(ax, 'nextplot', 'add');
plot(ax, Yeucl(:,2), 'color', cols{2}, 'linewidth', lw)
plot(ax, Yeucl(:,3), 'color', cols{3}, 'linewidth', lw)
plot(ax, [1 nS], [0 0], ':k', 'linewidth', lw)
ana.botmpaper.axesConfig1(ax, nC);

ax = axs(8);
plot(ax, Yconv(:,1), 'color', cols{1}, 'linewidth', lw)
set(ax, 'nextplot', 'add');
plot(ax, Yconv(:,2), 'color', cols{2}, 'linewidth', lw)
plot(ax, Yconv(:,3), 'color', cols{3}, 'linewidth', lw)
plot(ax, [1 nS], [0 0], ':k', 'linewidth', lw)
ana.botmpaper.axesConfig1(ax, nC);

ax = axs(9);
plot(ax, Ymatch(:,1), 'color', cols{1}, 'linewidth', lw)
set(ax, 'nextplot', 'add');
plot(ax, Ymatch(:,2), 'color', cols{2}, 'linewidth', lw)
plot(ax, Ymatch(:,3), 'color', cols{3}, 'linewidth', lw)
plot(ax, [1 nS], [0 0], ':k', 'linewidth', lw)
ana.botmpaper.axesConfig1(ax, nC);

ax = axs(10);
plot(ax, Ybotm(:,1), 'color', cols{1}, 'linewidth', lw)
set(ax, 'nextplot', 'add');
plot(ax, Ybotm(:,2), 'color', cols{2}, 'linewidth', lw)
plot(ax, Ybotm(:,3), 'color', cols{3}, 'linewidth', lw)
plot(ax, [1 nS], [0 0], ':k', 'linewidth', lw)
ana.botmpaper.axesConfig1(ax, nC);

linkaxes(axs(1:5), 'x');
set(axs(1:5), 'xlim', [1900    2850])
linkaxes(axs(6:10), 'x');
set(axs(6:10), 'xlim', [15448    16678])
% mysort.plot.savefig(f, 'B', 'ai', 1, 'fig', 0);

%% NEW PANEL FOR REVIEWER
%%
Yeucl  = zeros(nS, nT);
Yconv  = zeros(nS, nT);
Ymatch = zeros(nS, nT);
Ybotm = zeros(nS, nT);
for t=1:3
    T = squeeze(GT.XI(:,1,t));
    Yconv(:,t)  = conv2(GT.X(1:nS)', flipud(T)/norm(T), 'same');
    CC = (GT.C + GT.C(1,1)*eye(Tf))/2;
    F = (T'/CC)';
%     F = (T'*inv(GT.C))';
    Ymatch(:,t) = conv2(GT.X(1:nS)', flipud(F), 'same');
    Ybotm(:,t)  = Ymatch(:,t) - .5*F'*T + log(.0001);
    Ymatch(:,t) = Ymatch(:,t)/norm(F);
    Tf2 = round(Tf/2);
    for s=Tf2:nS-Tf2
        s1 = s-Tf2+1;
        s2 = s1+Tf-1;
        Yeucl(s,t) = norm(GT.X(s1:s2)'-T)/nC;
    end
end
% PLOT
f = open('B.fig');

axs = mysort.plot.getAxesWithTagFromFigure(f, 'axes1');
axs(2) = mysort.plot.getAxesWithTagFromFigure(f, 'axes6');
axs(3) = mysort.plot.getAxesWithTagFromFigure(f, 'axes7');
axs(4) = mysort.plot.getAxesWithTagFromFigure(f, 'axes8');
axs(5) = mysort.plot.getAxesWithTagFromFigure(f, 'axes9');
axs(6) = mysort.plot.getAxesWithTagFromFigure(f, 'axes10');
axs(7) = mysort.plot.getAxesWithTagFromFigure(f, 'axes11');
axs(8) = mysort.plot.getAxesWithTagFromFigure(f, 'axes12');
axs(9) = mysort.plot.getAxesWithTagFromFigure(f, 'axes13');
axs(10) = mysort.plot.getAxesWithTagFromFigure(f, 'axes14');

plotIDs = [4 20 8];
pgdf = gdf;
pgdf(:,1) = plotIDs(gdf(:,1));
pgdf2 = gdf2;
pgdf2(:,1) = plotIDs(gdf2(:,1));

cols = {mysort.plot.vectorColor(plotIDs(1)) mysort.plot.vectorColor(plotIDs(2)) mysort.plot.vectorColor(plotIDs(3))};
lw = 2;


ax = axs(4);
plot(ax, Ybotm(:,1), 'color', cols{1}, 'linewidth', lw)
set(ax, 'nextplot', 'add');
plot(ax, Ybotm(:,2), 'color', cols{2}, 'linewidth', lw)
plot(ax, Ybotm(:,3), 'color', cols{3}, 'linewidth', lw)
plot(ax, [1 nS], [0 0], ':k', 'linewidth', lw)
ana.botmpaper.axesConfig1(ax, nC);

[mx, mxidx] = max(Ybotm, [], 2);
Yamax = nan(size(Ybotm));
Yamax(mxidx==1,1) = Ybotm(mxidx==1,1);
Yamax(mxidx==2,2) = Ybotm(mxidx==2,2);
Yamax(mxidx==3,3) = Ybotm(mxidx==3,3);

ax = axs(5);
plot(ax, Yamax(:,1), 'color', cols{1}, 'linewidth', lw)
set(ax, 'nextplot', 'add');
plot(ax, Yamax(:,2), 'color', cols{2}, 'linewidth', lw)
plot(ax, Yamax(:,3), 'color', cols{3}, 'linewidth', lw)
plot(ax, [1 nS], [0 0], ':k', 'linewidth', lw)
ana.botmpaper.axesConfig1(ax, nC);


Yeucl  = zeros(nS, nT);
Yconv  = zeros(nS, nT);
Ymatch = zeros(nS, nT);
Ybotm = zeros(nS, nT);
for t=1:3
    T = squeeze(GT2.XI(:,1,t));
    Yconv(:,t)  = conv2(GT2.X(1:nS)', flipud(T)/norm(T), 'same');
    CC = (GT2.C + GT2.C(1,1)*eye(Tf))/2;
    F = (T'/CC)';
%     F = (T'*inv(GT.C))';
    Ymatch(:,t) = conv2(GT2.X(1:nS)', flipud(F), 'same');
    Ybotm(:,t)  = Ymatch(:,t) - .5*F'*T + log(.0001);
    Ymatch(:,t) = Ymatch(:,t)/norm(F);
    Tf2 = round(Tf/2);
    for s=Tf2:nS-Tf2
        s1 = s-Tf2+1;
        s2 = s1+Tf-1;
        Yeucl(s,t) = norm(GT2.X(s1:s2)'-T)/nC;
    end
end

ax = axs(9);
plot(ax, Ybotm(:,1), 'color', cols{1}, 'linewidth', lw)
set(ax, 'nextplot', 'add');
plot(ax, Ybotm(:,2), 'color', cols{2}, 'linewidth', lw)
plot(ax, Ybotm(:,3), 'color', cols{3}, 'linewidth', lw)
plot(ax, [1 nS], [0 0], ':k', 'linewidth', lw)
ana.botmpaper.axesConfig1(ax, nC);

ax = axs(10);
[mx, mxidx] = max(Ybotm, [], 2);
Yamax = nan(size(Ybotm));
Yamax(mxidx==1,1) = Ybotm(mxidx==1,1);
Yamax(mxidx==2,2) = Ybotm(mxidx==2,2);
Yamax(mxidx==3,3) = Ybotm(mxidx==3,3);
plot(ax, Yamax(:,1), 'color', cols{1}, 'linewidth', lw)
set(ax, 'nextplot', 'add');
plot(ax, Yamax(:,2), 'color', cols{2}, 'linewidth', lw)
plot(ax, Yamax(:,3), 'color', cols{3}, 'linewidth', lw)
plot(ax, [1 nS], [0 0], ':k', 'linewidth', lw)
ana.botmpaper.axesConfig1(ax, nC);

linkaxes(axs(1:5), 'x');
set(axs(1:5), 'xlim', [1900    2850])
linkaxes(axs(6:10), 'x');
set(axs(6:10), 'xlim', [15448    16678])
% mysort.plot.savefig(f, 'FigureForReviewerArgMaxOverTime_', 'ai', 1, 'fig', 0);
% mysort.plot.savefig(f, 'C:\Users\frankef\Dropbox\PaperBOTM1\SubmissionJCompNeuro\Revision1\FigureForReviewerArgMaxOverTime_', 'ai', 1, 'fig', 0);