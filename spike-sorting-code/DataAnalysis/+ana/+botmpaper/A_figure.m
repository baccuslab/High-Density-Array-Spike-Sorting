% f = mysort.plot.figure('w', 1200, 'h', 800);
% mysort.plot.savefig(f, 'A', 'png', 0);
% guide('A.fig')

%T = ana.torbjorn.websiteSlamander.getGTTemplates();
%save('torbjoernTemplates.mat', 'T');
load('torbjoernTemplates.mat');
T = mysort.wf.tAlignOnMax(-T); T = -T;
T = mysort.wf.tResample(T, 1, 3);
% elPos = ana.torbjorn.websiteSlamander.getElPos();
% mysort.plot.templates2D(T*5, elPos(:,2:3), 5, [], 'stacked',0)
[mi ma mi_idx ma_idx] = mysort.wf.tMinMaxPerTemplate(T);
%%
ID1 = 1; ID2 = 2;
[mins minChans] = sort(mi(ID1,:));
T1 = T(:,minChans(1:4),ID1); T1= T1/max(abs(T1(:)));
% T1(:,2) = 2*T1(:,2);
[mins minChans] = sort(mi(ID2,:));
T2 = T(:,minChans([2 1 3 4]),ID2); T2= T2/max(abs(T2(:)));
% T2(:,1) = 2*T2(:,1);
nS = 10000;
nSshow = 500;
nC = 4;
Tf = size(T,1);
noiseStd = .2;
X = noiseStd*randn(nS, nC);

vT1 = mysort.wf.t2v(T1);
vT2 = mysort.wf.t2v(T2);

%% ARTEFACT
A1L = 100;
% cycles = 6;
cycles = 5;
% phaseShiftBetweenChannels = round((A1L/cycles)/4);
% phaseShiftBetweenChannels = round( (A1L/cycles)/2 );
phaseShiftBetweenChannels = 0;

A1 = sin(.5*2*pi*(0:A1L-1)/(A1L)).*sin(cycles*2*pi*(0:A1L-1)/(A1L));
% figure; plot(A1)
A1t = [300 1000:200:8000];
artefactAmplitudesPerChannel = [1 2 3 1]/3;
% artefactAmplitudesPerChannel = [2 2 2 2];
for i=1:length(A1t)
    for c=1:nC
        phaseShift = c * phaseShiftBetweenChannels;
        s1 = A1t(i)+phaseShift;
        s2 = s1+A1L-1;
        X(s1:s2,c) = X(s1:s2,c) + A1'*artefactAmplitudesPerChannel(c);
    end
end
%%

t1 = 40+[20 140];
t2 = 40+[60 150 325];
for i=1:length(t1)
    X(t1(i):t1(i)+Tf-1,:) = X(t1(i):t1(i)+Tf-1,:) + T1;
end
for i=1:length(t2)
    X(t2(i):t2(i)+Tf-1,:) = X(t2(i):t2(i)+Tf-1,:) + T2;
end

XCC = mysort.noise.XCorrContainer(X, Tf-1);
f1 = mysort.util.embedChan2embedTime(XCC.invMul(mysort.util.embedTime2embedChan(vT1,nC), [1:nC]),nC);
f2 = mysort.util.embedChan2embedTime(XCC.invMul(mysort.util.embedTime2embedChan(vT2,nC), [1:nC]),nC);
tf1 = mysort.wf.v2t(f1,nC);
tf2 = mysort.wf.v2t(f2,nC);

Yconv1 = zeros(size(X));
Yeucl1 = zeros(size(X));

Ymatc1 = zeros(size(X));
Yconv2 = zeros(size(X));
Yeucl2 = zeros(size(X));
Ymatc2 = zeros(size(X));

for c=1:nC
    Yconv1(:,c) = conv2(X(:,c), flipud(T1(:,c)/norm((T1(:)))), 'same');
    Yconv2(:,c) = conv2(X(:,c), flipud(T2(:,c)/norm((T2(:)))), 'same');
    Ymatc1(:,c) = conv2(X(:,c), flipud(tf1(:,c)/norm((tf1(:)))), 'same');
    Ymatc2(:,c) = conv2(X(:,c), flipud(tf2(:,c)/norm((tf2(:)))), 'same');
    for t=1:nS
        s1 = t;
        s2 = min(nS, t+Tf-1);
        Yeucl1(t,c) = norm(X(s1:s2,c)-T1(1:s2-s1+1,c))/nC;
        Yeucl2(t,c) = norm(X(s1:s2,c)-T2(1:s2-s1+1,c))/nC;
    end
end
%%
if 0
    Ymaha1 = zeros(size(X));
    Ymaha1total = zeros(size(X,1), 1);
    Ymaha2 = zeros(size(X));
    Ymaha2total = zeros(size(X,1), 1);
    for t=1:nS-Tf
        s1 = t;
        s2 = min(nS, t+Tf-1);    
        d1 = mysort.wf.t2v(X(s1:s2,:)-T1(1:s2-s1+1,:));
        d2 = mysort.wf.t2v(X(s1:s2,:)-T2(1:s2-s1+1,:));
        n1 = mysort.wf.v2t(mysort.util.embedChan2embedTime(XCC.invMul(mysort.util.embedTime2embedChan(d1,nC), [1:nC]),nC),nC);
        n2 = mysort.wf.v2t(mysort.util.embedChan2embedTime(XCC.invMul(mysort.util.embedTime2embedChan(d2,nC), [1:nC]),nC),nC);
        for c=1:nC
            Ymaha1(t,c) = norm(n1(:,c));
            Ymaha2(t,c) = norm(n2(:,c));
        end
        Ymaha1total(t) = norm(n1);
        Ymaha2total(t) = norm(n2);
    end

    figure
    plot(sum(Yeucl1(1:nSshow,:),2)')
    hold on
    plot(sum(Yeucl2(1:nSshow,:),2)', 'r')
    plot(sum(Ymaha1total(1:nSshow,:),2)'/7, 'b:')
    plot(sum(Ymaha2total(1:nSshow,:),2)'/7, 'r:')
    
    
end

if 1
    %%
    figure
    plot(sum(Yconv1,2))
    hold on
    plot( mysort.util.mcfilt(X', mysort.wf.v2m(f1,nC)/norm(f1),'same'), 'g');
end

vX = mysort.wf.t2v(X(1:nSshow,:));
vT1 = mysort.wf.t2v(T1);
vT2 = mysort.wf.t2v(T2);
vYconv1 = mysort.wf.t2v(Yconv1(1:nSshow,:));
vYeucl1 = mysort.wf.t2v(Yeucl1(1:nSshow,:));
vYmatc1 = mysort.wf.t2v(Ymatc1(1:nSshow,:));
vYconv2 = mysort.wf.t2v(Yconv2(1:nSshow,:));
vYeucl2 = mysort.wf.t2v(Yeucl2(1:nSshow,:));
vYmatc2 = mysort.wf.t2v(Ymatc2(1:nSshow,:));

f = open('A.fig');
cSpacer = 2;
%% Data
IDdata = 18;
IDT1 = 15;
IDT2 = 16;
ax = mysort.plot.getAxesWithTagFromFigure(f, 'axes1');
mysort.plot.waveformsVertical(vX, 'nC', nC, 'axesHandle', ax, 'channelSpacer', repmat(cSpacer,1,nC), 'IDs', IDdata)
ana.botmpaper.axesConfig1(ax, nC);

ax = mysort.plot.getAxesWithTagFromFigure(f, 'axes2');
mysort.plot.waveformsVertical(vX, 'nC', nC, 'axesHandle', ax, 'channelSpacer', repmat(cSpacer,1,nC), 'IDs', IDdata)
ana.botmpaper.axesConfig1(ax, nC);

ax = mysort.plot.getAxesWithTagFromFigure(f, 'axes11');
mysort.plot.waveformsVertical(vX, 'nC', nC, 'axesHandle', ax, 'channelSpacer', repmat(cSpacer,1,nC), 'IDs', IDdata)
ana.botmpaper.axesConfig1(ax, nC);

%% Templates
ax = mysort.plot.getAxesWithTagFromFigure(f, 'axes3');
mysort.plot.waveformsVertical([vT1; vT2], 'nC', nC, 'axesHandle', ax, 'channelSpacer', repmat(cSpacer,1,nC), 'IDs', [IDT1 IDT2], 'linewidth', 2)
ana.botmpaper.axesConfig1(ax, nC);

ax = mysort.plot.getAxesWithTagFromFigure(f, 'axes4');
mysort.plot.waveformsVertical([vT1; vT2], 'nC', nC, 'axesHandle', ax, 'channelSpacer', repmat(cSpacer,1,nC), 'IDs', [IDT1 IDT2], 'linewidth', 2)
ana.botmpaper.axesConfig1(ax, nC);

ax = mysort.plot.getAxesWithTagFromFigure(f, 'axes12');
mysort.plot.waveformsVertical([f1; f2]/max(max(abs([f1; f2]))), 'nC', nC, 'axesHandle', ax, 'channelSpacer', repmat(cSpacer,1,nC), 'IDs', [IDT1 IDT2], 'linewidth', 2)
ana.botmpaper.axesConfig1(ax, nC);

%% Conv
ax = mysort.plot.getAxesWithTagFromFigure(f, 'axes7');
mysort.plot.waveformsVertical([vYeucl1; vYeucl2]*1.5, 'nC', nC, 'axesHandle', ax, 'channelSpacer', repmat(cSpacer,1,nC), 'IDs', [IDT1 IDT2], 'linewidth', 2)
ana.botmpaper.axesConfig1(ax, nC);

ax = mysort.plot.getAxesWithTagFromFigure(f, 'axes8');
mysort.plot.waveformsVertical([vYconv1; vYconv2]/1.5, 'nC', nC, 'axesHandle', ax, 'channelSpacer', repmat(cSpacer,1,nC), 'IDs', [IDT1 IDT2], 'linewidth', 2)
ana.botmpaper.axesConfig1(ax, nC);

ax = mysort.plot.getAxesWithTagFromFigure(f, 'axes13');
mysort.plot.waveformsVertical([vYmatc1; vYmatc2]/1.5, 'nC', nC, 'axesHandle', ax, 'channelSpacer', repmat(cSpacer,1,nC), 'IDs', [IDT1 IDT2], 'linewidth', 2)
ana.botmpaper.axesConfig1(ax, nC);


%% Result
ax = mysort.plot.getAxesWithTagFromFigure(f, 'axes9');
mysort.plot.waveformsVertical(sum(Yeucl1(1:nSshow,:),2)', 'nC', 1, 'axesHandle', ax, 'IDs', IDT1, 'linewidth', 2)
mysort.plot.waveformsVertical(sum(Yeucl2(1:nSshow,:),2)', 'nC', 1, 'axesHandle', ax, 'IDs', IDT2, 'linewidth', 2)
ana.botmpaper.axesConfig1(ax, 1);

ax = mysort.plot.getAxesWithTagFromFigure(f, 'axes10');
mysort.plot.waveformsVertical(sum(Yconv1(1:nSshow,:),2)', 'nC', 1, 'axesHandle', ax, 'IDs', [IDT1], 'linewidth', 2)
hold on
mysort.plot.waveformsVertical(sum(Yconv2(1:nSshow,:),2)', 'nC', 1, 'axesHandle', ax, 'IDs', [IDT2], 'linewidth', 2)
ana.botmpaper.axesConfig1(ax, 1);

ax = mysort.plot.getAxesWithTagFromFigure(f, 'axes14');
mysort.plot.waveformsVertical(sum(Ymatc1(1:nSshow,:),2)', 'nC', 1, 'axesHandle', ax, 'IDs', IDT1, 'linewidth', 2)
mysort.plot.waveformsVertical(sum(Ymatc2(1:nSshow,:),2)', 'nC', 1, 'axesHandle', ax, 'IDs', IDT2, 'linewidth', 2)
ana.botmpaper.axesConfig1(ax, 1);

%% Export
% mysort.plot.savefig(f, 'A', 'fig', 0, 'ai', 1);