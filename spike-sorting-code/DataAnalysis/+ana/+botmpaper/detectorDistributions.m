load('torbjoernTemplates.mat');
T = mysort.wf.tAlignOnMax(-T); T = -T;
T = mysort.wf.tResample(T, 2, 3);
templatePeakHeight = 4;
[mi ma mi_idx ma_idx] = mysort.wf.tMinMaxPerTemplate(T);
ID1 = 1; ID2 = 2;
[mins minChans] = sort(mi(ID1,:));
T1 = T(:,minChans(1:4),ID1); T1 = templatePeakHeight*T1/max(abs(T1(:)));
% T1(:,2) = 2*T1(:,2);
[mins minChans] = sort(mi(ID2,:));
T2 = T(:,minChans([2 1 3 4]),ID2); T2 = templatePeakHeight*T2/max(abs(T2(:)));
vT1 = mysort.wf.t2v(T1);
vT2 = mysort.wf.t2v(T2);
[mi ma mi_idx ma_idx] = mysort.wf.tMinMaxPerTemplate(T1);
peakPositionInTemplate1 = mi_idx(1);
[mi ma mi_idx ma_idx] = mysort.wf.tMinMaxPerTemplate(T2);
peakPositionInTemplate2 = mi_idx(1);

figure;
plot([vT1' vT2'])
hold on
plot(peakPositionInTemplate1, vT1(peakPositionInTemplate1), 'rd');
plot(peakPositionInTemplate2, vT2(peakPositionInTemplate1), 'rd');

Tf = size(T,1);
nC = 4; 
noiseStd = 1;

nS = 100000;
X = noiseStd*randn(nS, nC*Tf);

% Eucl:
eucl_sig_plus_noise = sum(X.^2, 2);
eucl_noise_only = sum((X-repmat(vT1,nS,1)).^2, 2);

% conv:
conv_sig_plus_noise = X*vT1' + vT1*vT1';
conv_noise_only = X*vT1';

%
m_eucl = min(min(eucl_sig_plus_noise),min(eucl_noise_only));
m_conv = min(min(conv_sig_plus_noise),min(conv_noise_only));
M_eucl = max(max(eucl_sig_plus_noise),max(eucl_noise_only));
M_conv = max(max(conv_sig_plus_noise),max(conv_noise_only));

edges = linspace(min(m_eucl, m_conv), max(M_eucl, M_conv), 200);
binWidth = edges(2)-edges(1);
n=[];
n(:,1) = histc(eucl_sig_plus_noise, edges)/(nS*binWidth);
n(:,2) = histc(eucl_noise_only, edges)/(nS*binWidth);
n(:,3) = histc(conv_sig_plus_noise, edges)/(nS*binWidth);
n(:,4) = histc(conv_noise_only, edges)/(nS*binWidth);
% n(n==0) = nan;

%%
dims = Tf*nC;
% dims = 5;
% for a given number of dimensions compute the energy of the signal only on
% so many samples as dimensions are available:
[sT1 t1 t2] = ana.botmpaper.cutTemplateToMiddleWithNDims(T1, peakPositionInTemplate2, Tf);
e = norm(sT1);

scalefactor = 1/noiseStd;
ncx = scalefactor^2*ncx2pdf(scalefactor^2*edges, dims, (scalefactor*e)^2);
ccx = scalefactor^2*chi2pdf(scalefactor^2*edges, dims);

% Compute overlap of distributions
NoiseDist  = @(x) ncx2cdf(scalefactor^2*x, dims, (scalefactor*e)^2);
SignalDist = @(x) chi2cdf(scalefactor^2*x, dims);
pNoise = 10000;
ErrorF = @(x) (pNoise*NoiseDist(x) - SignalDist(x));
opt.Display = 'on';
[xstar fxstar] = fminsearch(ErrorF, dims/scalefactor^2);
figure
xr = .5*dims/scalefactor^2:.1:1.5*dims/scalefactor^2;
plot(xr, NoiseDist(xr), 'r');
hold on
plot(xr, SignalDist(xr), 'g');
plot(xr, ErrorF(xr), 'b');
plot(xstar, fxstar, 'md');
set(gca, 'ylim', [-1 1])
title(sprintf('dims = %d', dims))

%%

bestVals1 = [];
bestXVals1 = [];
bestVals2 = [];
bestXVals2 = [];
TfRange = [1:1:40];
xr = -10:.01:50;
% figure
% ah = subplot(2,1,1);
% ah(2) = subplot(2,1,2);
% scalefactor = 1;
for tfi=1:length(TfRange)
    Tf = TfRange(tfi);
    d = Tf*nC;
    
    
    NoiseDist  = @(x, xe) ncx2cdf(scalefactor^2*x, d, (scalefactor*xe)^2);
    SignalDist = @(x) chi2cdf(scalefactor^2*x, d);
    pNoise = 1000;
    ErrorF = @(x, xe) (pNoise*NoiseDist(x, xe) - SignalDist(x));
    [sT1 t1 t2] = ana.botmpaper.cutTemplateToMiddleWithNDims(T1, peakPositionInTemplate1, Tf);
    [sT2 t1 t2] = ana.botmpaper.cutTemplateToMiddleWithNDims(T2, peakPositionInTemplate2, Tf);
    e1 = norm(sT1);
    e2 = norm(sT2);    

    x0 = d/scalefactor^2;
%     x0 = d;
    [bestXVals1(tfi) bestVals1(tfi)] = fminsearch(@(x) ErrorF(x,e1), x0);
    [bestXVals2(tfi) bestVals2(tfi)] = fminsearch(@(x) ErrorF(x,e2), x0);
%     cla(ah(1))
%     plot(ah(1), xr, NoiseDist(xr, e1), 'r');
%     set(ah(1), 'nextplot', 'add');
%     plot(ah(1), xr, SignalDist(xr), 'g');
%     plot(ah(1), x0, 0, 'om');
%     
%     cla(ah(2))
%     plot(ah(2), xr, ErrorF(xr,e1));
%     set(ah(2), 'nextplot', 'add');
%     plot(ah(2), bestXVals1(tfi), 0, 'rd');
%     pause
end
%%
figure;
plot(ErrorF(0:.1:100,e1));
%%
figure;
subplot(2,1,1)
plot(TfRange,bestVals1)
subplot(2,1,2)
plot(TfRange,bestXVals1)

%%
fh = mysort.plot.figure('w', 900, 'h', 450);
ah = axes;
set(ah, 'fontsize', 14);
plot(TfRange, -100*bestVals1, 'linewidth', 2, 'color', mysort.plot.vectorColor(15))
hold on
plot(TfRange, -100*bestVals2, 'linewidth', 2, 'color', mysort.plot.vectorColor(16))
xlabel('Template dimensions (nC*L) [samples]')
ylabel('Detection Performance [%]')
legend('Template 1', 'Template 2')

%%
fh = mysort.plot.figure('w', 900, 'h', 450);
ah = subplot(1,3,1);
set(ah, 'fontsize', 14);

noisecolor = [.8 0 0];
signalcolor = [0 .8 0];

%h = bar(ah(1), edges, n(:,[2 1]), 'stacked');
% set(h(1), 'facecolor', noisecolor, 'edgecolor', noisecolor)
% set(h(2), 'facecolor', signalcolor, 'edgecolor', signalcolor)
h = plot(ah(1), edges+binWidth/2, n(:,[2 1]), 'linewidth', 2);
hold on
% plot(ah(1), edges, ncx, ':b', 'linewidth', 2);
% plot(ah(1), edges, ccx, ':g', 'linewidth', 2);


set(h(1), 'color', noisecolor)
set(h(2), 'color', signalcolor)
xlabel('detector output');
ylabel('pdf');
title('euclidean distance')


ah(2) = subplot(1,3,2);
set(ah(2), 'fontsize', 14);

h = plot(ah(2), edges+binWidth/2, n(:,[4 3]), 'linewidth', 2);
set(h(1), 'color', noisecolor)
set(h(2), 'color', signalcolor)
xlabel('detector output');

% ylabel('pdf');
legend('noise only', 'template + noise')
title('convolution detector')

set(ah(1), 'xlim', [100 450]); %'ylim', [0 .15],
set(ah(2), 'ylim', [0 .05], 'xlim', [-50 150]);

ah(2) = subplot(1,3,2);
set(ah(2), 'fontsize', 14);


% mysort.plot.savefig(fh, 'detectorDistributions', 'png', 0, 'fig', 0, 'ai', 1)
