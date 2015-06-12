close all
T = [];
L = 5;
nE = 3;
T(:,:,1) = [0 0 0
            1 2 0
            0 0 0];
T(:,:,2) = [0 0 0
            1 2 1
            0 0 0];
T(:,:,3) = [0 0 0
            2 2 0
            0 0 0];
T = T*3;
priors = [.01 .02 .03];


elSet = [1 2];
[Q0, E, thr, P] = mysort.wf.qualityOfTemplateSet(T, priors, true);
mysort.plot.figureTitle(sprintf('Els: all, Q: %.3f', Q0));

elSet = [1 2];
[Q1, E, thr, P] = mysort.wf.qualityOfTemplateSet(T(:,elSet,:), priors, true);
mysort.plot.figureTitle(sprintf('Els: %d %d, Q: %.3f', elSet, Q1));

elSet = [1 3];
[Q2, E, thr, P] = mysort.wf.qualityOfTemplateSet(T(:,elSet,:), priors, true);
mysort.plot.figureTitle(sprintf('Els: %d %d, Q: %.3f', elSet, Q2));

elSet = [2 3];
[Q3, E, thr, P] = mysort.wf.qualityOfTemplateSet(T(:,elSet,:), priors, true);
mysort.plot.figureTitle(sprintf('Els: %d %d, Q: %.3f', elSet, Q3));