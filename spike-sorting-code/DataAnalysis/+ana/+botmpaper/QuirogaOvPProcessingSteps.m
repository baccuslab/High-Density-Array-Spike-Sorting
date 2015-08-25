benchmarks = {'Easy1', 'Easy2', 'Difficult1', 'Difficult2'};
noiselevels = {'005', '01', '015', '02'};
initTf = 81;
initCutLeft = -5;
Tf = 61;
cutLeft = 15;
dpath = ana.botmpaper.E10_quirogaDataPath();
samplesPerSecond = 32000;
Eucl  = @ana.botmpaper.Eucl;
Maha  = @ana.botmpaper.Maha;
Conv  = @ana.botmpaper.Conv;
Match = @ana.botmpaper.Match;
Botm  = @ana.botmpaper.Botm;
global savefig__;
savefig__ = 1;
figSavePath = 'C:\Dropbox\PaperBOTM1\DetectorFigures\';

%%
RR = [];
b = 4 %length(benchmarks)
n = 2 %length(noiselevels)
savePath = fullfile(dpath, benchmarks{b}, noiselevels{n});
if ~exist(savePath, 'file')
    mkdir(savePath);
end      

pic = 0;
fprintf('Starting with %s %s\n', benchmarks{b}, noiselevels{n});
close all
quirogaFile = sprintf('C_%s_noise%s.mat', benchmarks{b}, noiselevels{n});        
GT = ana.botmpaper.E10_preprocessing(dpath, quirogaFile, initCutLeft, initTf, dpath);

bufferFile = fullfile(savePath, 'botmOutputDistributionsThrehold.mat');
        
% RUN QuirogaOptimalThreshold first to save all buffered files.
load(bufferFile);
%% COMPUTE FILTER OUTPUTS
Y = zeros(size(GT.X));
Tf2 = floor(size(T,2)/2);
M = zeros(2*Tf2+1, 3, 3);
for i=1:3
    t = T(i,:);
    f = t/C_;
    Y(i,:)  = conv(GT.X, flipud(f(:)), 'same') - .5*t*f';
    for j=1:3
        M(:,i,j) = conv(T(j,:), flipud(f(:)), 'same');
    end
end
maxTau = 4;
DC = mysort.sorters.DiscriminantFunctionContainer(Y', M, maxTau);

%% annotation of peaks
OVERLAP_DIST = 15;
OTHERSPIKE_DIST = 6;
SINGLESPIKE_DIST = 3;
D = [];
for i=1:3
    p = peaks{i};
    p(end, 3) = 0;
    mySpikeTimes = tgdf(tgdf(:,1)==i,2) + cutLeft + 1;
    otherSpikeTimes = tgdf(tgdf(:,1)~=i,2) + cutLeft + 1;
    for k=1:size(p,1)
        dists = abs(mySpikeTimes - p(k,1));
        di = find(dists < SINGLESPIKE_DIST);
        do = abs(otherSpikeTimes - p(k,1));
        bIsOvp = any(do < OVERLAP_DIST);
        bIsOtherSpike = any(do < OTHERSPIKE_DIST);
        if ~isempty(di);
            D(end+1) = dists(di);
            if bIsOvp
                p(k,3) = 2; % ovp
            else
                p(k,3) = 1; % single spike
            end
        else
            if bIsOtherSpike
                p(k,3) = 3; % other spike
            else
                p(k,3) = 4; % noise
            end
        end
    end
    peaks{i} = p;
end

%% SELECT AN OVERLAPPING SPIKE
ovpIdx = find(p(:,3) == 2);
ovpIdx = ovpIdx(2);
% find spikes in gt
s_idx = find(abs(tgdf(:,2)-p(ovpIdx,1)) < 25);
% s_idx = [s_idx s_idx(1)+1];
s1_id = tgdf(s_idx(1),1);
s1_t  = tgdf(s_idx(1),2);
s2_id = tgdf(s_idx(2),1);
s2_t  = tgdf(s_idx(2),2);
dt = s2_t-s1_t;
region = 160;
L = size(T,2) + dt + region;

[maxPerOvpFun maxTimeIdxPerOvpFun] = max(DC.DOvp(start:start+L,:));
[sortedMaxima maximaIdx] = sort(maxPerOvpFun, 'descend');
DC.DOvpIndex(maximaIdx(1:3),:)
cols = {'b', 'r', 'g'};

% Compute Error over threhold
fh2 = mysort.plot.figure('w',300, 'h', 900);
ph=0; ah=[]; nR= 8; nC= 1;
for i=1:3
    ph=ph+1; ah(ph) = subplot(nR, nC, ph);
    plot(T(i,:), cols{i}, 'linewidth', 2)
    set(ah(ph), 'xlim', [0 L]-region/2, 'ylim', [-1 2]);
end

ph=ph+1; ah(ph) = subplot(nR, nC, ph);
start = s1_t-cutLeft-region/2;
plot(GT.X(start:start+L), 'k', 'linewidth', 2)
set(ah(ph), 'xlim', [0 L], 'ylim', [-1 2]);

ph=ph+1; ah(ph) = subplot(nR, nC, ph); 
start = s1_t-cutLeft-region/2;
hold on
for i=1:3
    plot(Y(i,start:start+L)', cols{i}, 'linewidth', 2)
end
plot([0 L], [0 0], 'k:', 'linewidth', 2)
set(ah(ph), 'xlim', [0 L]);

ph=ph+1; ah(ph) = subplot(nR, nC, ph); 
hold on
plot(DC.DOvp(start:start+L,:), 'color', [.7 .7 .85], 'linewidth', 1)
for i=1:3
    plot(Y(i,start:start+L)', cols{i}, 'linewidth', 2)
end
plot(DC.DOvp(start:start+L,maximaIdx(1)), 'color', [.2 .2 .4], 'linewidth', 2)

plot([0 L], [0 0], 'k:', 'linewidth', 2)
set(ah(ph), 'xlim', [0 L]);

SIC0 = resample(Y(:,start:start+L)', 3, 1)';
ph=ph+1; ah(ph) = subplot(nR, nC, ph); 
hold on
for i=1:3
    plot(SIC0(i,:), cols{i}, 'linewidth', 2)
end
L = L*3;
[mx maxTemplate] = max(SIC0);
[mx mxidx] = max(mx);
maxTemplate = maxTemplate(mxidx);
plot(mxidx, 1.1*mx, '*', 'color', cols{maxTemplate}, 'linewidth', 1)
plot([0 L], [0 0], 'k:', 'linewidth', 2)
set(ah(ph), 'xlim', [0 L]);

SIC1 = SIC0;
Subtractor = resample(squeeze(M(:,:,maxTemplate)), 3, 1)';
SIC1 = mysort.util.shiftSubtract(SIC1, Subtractor, mxidx-3*(Tf2)-1); 
ph=ph+1; ah(ph) = subplot(nR, nC, ph); 
hold on
for i=1:3
    plot(SIC1(i,:)', cols{i}, 'linewidth', 2)
end
[mx maxTemplate] = max(SIC1);
[mx mxidx] = max(mx);
maxTemplate = maxTemplate(mxidx);
plot(mxidx, 1.1*mx, '*', 'color', cols{maxTemplate}, 'linewidth', 1)
plot([0 L], [0 0], 'k:', 'linewidth', 2)
set(ah(ph), 'xlim', [0 L], 'ylim', get(ah(ph-1), 'ylim'));

set(ah, 'xticklabel', [], 'yticklabel', [], 'xcolor', 'w', 'ycolor', 'w', 'box', 'off')

filePrefix = sprintf('%s_%s', benchmarks{b}, noiselevels{n});
mysort.plot.savefig(fh2, fullfile(figSavePath, [filePrefix 'flowchart_data']), 'ai', 1, 'fig', 0);

fprintf('Done.\n');
