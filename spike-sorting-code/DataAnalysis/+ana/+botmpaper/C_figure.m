%% Prepare Figure (only once)
% f = mysort.plot.figure('w', 1200, 'h', 800);
% mysort.plot.savefig(f, 'C', 'png', 0);
% guide('C.fig')

%% INIT

% %% Load Harris Eval
H = struct();
[H.REJECTIONPERF H.CLASSPERF H.DETPERF H.CPUTIMEPERF H.correctSpikesIdx H.otherSpikesIdx H.L] = ana.botmpaper.loadHarrisEval();

%% Load Quiroga Eval
Q = struct();
[Q.CLASSPERF Q.CPUTIMEPERF Q.DETPERF Q.L] = ana.botmpaper.loadQuirogaEval();

%% Restrict what to Plot
methodIdx = [1:3 5 6];
%                 1  2    3     4   5    6      7     8       9     10    11
preprocIdx = 11; %C, dC, C_opt, C_, CP, Cest, dCest, Cest_, CPest, CDL, CestDL
preprocIdx = 4; %C, dC, C_opt, C_, CP, Cest, dCest, Cest_, CPest, CDL
methodNames = Q.L.methods(methodIdx,4);
methodNames{end} = 'BOTM';
perfColor = [.0 .5 .1];

% Compute Quantities
QTotal = (Q.CLASSPERF + Q.DETPERF)/2;
QTotalM = squeeze(mean(QTotal));
QTotalStd = squeeze(std(QTotal));
QIndividualM = [squeeze(mean(Q.DETPERF(:,preprocIdx, methodIdx))) squeeze(mean(Q.CLASSPERF(:,preprocIdx, methodIdx)))];
QIndividualStd = [squeeze(std(Q.DETPERF(:,preprocIdx, methodIdx))) squeeze(std(Q.CLASSPERF(:,preprocIdx, methodIdx)))];
QCPUM = squeeze(mean(Q.CPUTIMEPERF(:,preprocIdx, methodIdx)));
QCPUStd = squeeze(std(Q.CPUTIMEPERF(:,preprocIdx, methodIdx)));

HTotal = (H.CLASSPERF + H.DETPERF + H.REJECTIONPERF)/3;
HTotalM = squeeze(mean(HTotal));
HTotalStd = squeeze(std(HTotal));
HIndividualM = [squeeze(mean(H.DETPERF(:,preprocIdx, methodIdx))) squeeze(mean(H.REJECTIONPERF(:,preprocIdx, methodIdx))) squeeze(mean(H.CLASSPERF(:,preprocIdx, methodIdx)))];
HIndividualStd = [squeeze(std(H.DETPERF(:,preprocIdx, methodIdx))) squeeze(std(H.REJECTIONPERF(:,preprocIdx, methodIdx))) squeeze(std(H.CLASSPERF(:,preprocIdx, methodIdx)))];
HCPUM = squeeze(mean(H.CPUTIMEPERF(:,preprocIdx, methodIdx)));
HCPUStd = squeeze(std(H.CPUTIMEPERF(:,preprocIdx, methodIdx)));

disp('Total Perf:')
[QTotalM(preprocIdx, methodIdx); HTotalM(preprocIdx, methodIdx)]

squeeze(H.CLASSPERF(:,4,:))
squeeze(H.CLASSPERF(:,preprocIdx,:))
%% PRINT CLASSIFICATION PERFORMANCES
% get average classification performance
classperfMQ = mean(squeeze(Q.CLASSPERF(:,preprocIdx,methodIdx)));
classperfMH = mean(squeeze(H.CLASSPERF(:,preprocIdx,methodIdx)));
% get error by euclidean template matching
errorEDQ = 100-classperfMQ(1);
errorEDH = 100-classperfMH(1);
% get error by MD matching
errorMDQ = 100-classperfMQ(2);
errorMDH = 100-classperfMH(2);
% compute performance increase
fprintf('Classification Error Q, ED: %f  MD: %f, %%decrease: %f \n', errorEDQ, errorMDQ, 100*(errorEDQ-errorMDQ)/errorEDQ);
fprintf('Classification Error H, ED: %f  MD: %f, %%decrease: %f \n', errorEDH, errorMDH, 100*(errorEDH-errorMDH)/errorEDH);

%%
assert(~any(size(H.L.methods)~=size(Q.L.methods)), 'Eval methods must be identical!')
f = open('C.fig');
axs = mysort.plot.getAxesWithTagFromFigure(f, 'axes1');
axs(2) = mysort.plot.getAxesWithTagFromFigure(f, 'axes2');
axs(3) = mysort.plot.getAxesWithTagFromFigure(f, 'axes3');
axs(4) = mysort.plot.getAxesWithTagFromFigure(f, 'axes4');
axs(5) = mysort.plot.getAxesWithTagFromFigure(f, 'axes5');
axs(6) = mysort.plot.getAxesWithTagFromFigure(f, 'axes6');
set(axs, 'fontsize', 14);

%--------------------------------------------------------------------------
% Total Performance Plots
bar(axs(1), QTotalM(preprocIdx, methodIdx), 'facecolor', perfColor)
set(axs(1), 'nextplot', 'add');
errorbar(axs(1), QTotalM(preprocIdx, methodIdx), QTotalStd(preprocIdx, methodIdx), 'k.', 'linewidth', 2, 'marker', 'none')
axes(axs(1));
matlabfilecentral.xticklabel_rotate.xticklabel_rotate(1:5,45,methodNames,'interpreter','none');
ylabel('Total Performance [%]')

bar(axs(2), HTotalM(preprocIdx, methodIdx), 'facecolor', perfColor)
set(axs(2), 'nextplot', 'add');
errorbar(axs(2), HTotalM(preprocIdx, methodIdx), HTotalStd(preprocIdx, methodIdx), 'k.', 'linewidth', 2, 'marker', 'none')
linkaxes(axs(1:2), 'y');
set(axs(2), 'yticklabel', []);
axes(axs(2));
matlabfilecentral.xticklabel_rotate.xticklabel_rotate(1:5,45,methodNames,'interpreter','none');

%--------------------------------------------------------------------------
% Individual Performance Plots
h1 = bar(axs(3), QIndividualM);
set(axs(3), 'nextplot', 'add');
for i=1:2
    if i==1
        off = -.15;
    else
        off = .15;
    end
    errorbar(axs(3), (1:5)+off, QIndividualM(:,i), QIndividualStd(:,i), 'k.', 'linewidth', 2, 'marker', 'none')
end
set(h1(1), 'faceColor', [0 0 1]);
set(h1(2), 'faceColor', [0 1 0]);
axes(axs(3));
matlabfilecentral.xticklabel_rotate.xticklabel_rotate(1:5,45,methodNames,'interpreter','none');
ylabel('Performance [%]')


h2 = bar(axs(4), HIndividualM);
set(axs(4), 'nextplot', 'add');
for i=1:3
    if i==1
        off = -.22;
    elseif i==2
        off = 0;
    else
        off = .22;
    end
    errorbar(axs(4), (1:5)+off, HIndividualM(:,i), HIndividualStd(:,i), 'k.', 'linewidth', 2, 'marker', 'none')
end
set(h2(1), 'faceColor', [0 0 1]);
set(h2(2), 'faceColor', [1 0 0]);
set(h2(3), 'faceColor', [0 1 0]);
linkaxes(axs(3:4), 'y');
set(axs(4), 'yticklabel', []);
axes(axs(4));
matlabfilecentral.xticklabel_rotate.xticklabel_rotate(1:5,45,methodNames,'interpreter','none');
h = legend('Detection', 'Rejection', 'Classification', 'location', 'east');
set(h,  'fontsize', 12); %'box', 'off',
set(axs(1:4), 'ylim', [0 110])

%--------------------------------------------------------------------------
% CPU Time Plots
bar(axs(5), QCPUM, 'facecolor', [.6 .6 .6])
set(axs(5), 'nextplot', 'add');
errorbar(axs(5), QCPUM, QCPUStd, 'k.', 'linewidth', 2, 'marker', 'none')
axes(axs(5));
matlabfilecentral.xticklabel_rotate.xticklabel_rotate(1:5,45,methodNames,'interpreter','none');
ylabel(axs(5), 'CPU Time [s]')
set(axs(5), 'ylim', [0 .013]);

bar(axs(6), HCPUM, 'facecolor', [.6 .6 .6])
set(axs(6), 'nextplot', 'add');
errorbar(axs(6), HCPUM, HCPUStd, 'k.', 'linewidth', 2, 'marker', 'none')
% linkaxes(axs(5:6), 'y');
% set(axs(6), 'yticklabel', []);
ylabel('CPU Time [s]')
axes(axs(6));
matlabfilecentral.xticklabel_rotate.xticklabel_rotate(1:5,45,methodNames,'interpreter','none');
set(axs(6), 'ylim', [0 .3]);



%% 
% mysort.plot.savefig(gcf, 'C', 'ai', 1, 'fig', 0);