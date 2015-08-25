mysort.plot.figure([800 300]);
%%
load('harrisMatchingNoisyTemplatesBuffer2.mat')

subplot(1,3,2);
% B = [100*std(Sensitivity,[],2) min(100, 100*(1-mean(Sensitivity,2)) + 100*std(Sensitivity,[],2))-100*(1-mean(Sensitivity,2))];
% boundedline(nSpikesPerTemplate, 100*mean(1-Sensitivity,2), B, 'g')
hold on

B = [100*std(Classification,[],2) 100*std(Classification,[],2)];
boundedline(nSpikesPerTemplate, 100*(mean(Classification,2)), B, 'b')

B = [100*std(Specificity,[],2) 100*std(Specificity,[],2)];
boundedline(nSpikesPerTemplate, 100*(1-mean(Specificity,2)), B, 'r')

box off
axis(gca, 'tight');
% set(gca, 'ylim', [30 100])
xlabel('% Of assigned Spikes Per Template')
ylabel('Performance [%]');
% hl = legend('Sensitivity', '', 'Classification',  'location', 'southwest'); set(hl, 'box', 'off');

ah = subplot(1,3,3);
P_ = mysort.plot.waveformsVertical(mysort.wf.v2t(T1([correctIdx],:), nC), 'axesHandle', ah, 'plotArgs', {'r', 'linewidth', 2}, 'ignoreClassesForColoring', 1, 'channelSpacer', 400); 
mysort.plot.waveformsVertical(mysort.wf.v2t(T2([correctIdx],:), nC), 'axesHandle', ah, 'plotArgs', {'b', 'linewidth', 2}, 'ignoreClassesForColoring', 1, 'channelSpacer', 400); 
set(ah, 'box', 'off')

%%
load('harrisMatchingDeletedTemplatesBuffer2.mat')

ah = subplot(1,3,1);
% B = [100*std(Sensitivity,[],2) 100*std(Sensitivity,[],2)];
% boundedline(nDeleteTemplates, 100*(1-mean(Sensitivity,2)), B, 'g')
hold on

B = [100*std(Classification,[],2) 100*std(Classification,[],2)];
plotMax = B(:,2) + 100*(mean(Classification,2));
plotMax = max(0, plotMax-100);
B(:,2) = B(:,2) - plotMax;
boundedline(nDeleteTemplates, 100*(mean(Classification,2)), B, 'b')

B = [100*std(Specificity,[],2) 100*std(Specificity,[],2)];
plotMax = B(:,2) + 100*(1-mean(Specificity,2));
plotMax = max(0, plotMax-100);
B(:,2) = B(:,2) - plotMax;
boundedline(nDeleteTemplates, 100*(1-mean(Specificity,2)), B, 'r')
box off
axis(gca, 'tight');
% set(gca, 'ylim', [70 100])
xlabel('#Deleted Templates')
ylabel('Performance [%]');
% hl = legend('Sensitivity', '', 'Classification', '', 'Specificity', 'location', 'southwest'); set(hl, 'box', 'off');
hl = legend('Classification', '', 'Specificity', 'location', 'southwest'); set(hl, 'box', 'off');
% mysort.plot.savefig(gcf, 'harrisMatchingNewFigure7', 'fig', 0, 'ai', 1);
% mysort.plot.savefig(gcf, 'C:\Users\frankef\Dropbox\PaperBOTM1\SubmissionJCompNeuro\Revision1\harrisMatchingNewFigure7', 'fig', 0, 'ai', 1);