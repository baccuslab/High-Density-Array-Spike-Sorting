sps = 20000;
name = 'TestMS';
nSessions = 5;
X = cell(1, nSessions);
sessionLength = [1000 2000 1300 4000 1500];
nChannels = 4;
for i=1:5
    X{i} = randn(sessionLength(i), nChannels)+2*i;
end
MSM = mysort.ds.MultiSessionMatrix(name, X, sps);

% mysort.plot.SliderDataAxes(MSM)
MSM.setActiveSession(4);

MSM2 = MSM.copy();
MSM2.restrictToChannels(1:2);
mysort.plot.SliderDataAxes({MSM MSM2})

size(MSM2)

MSM.concatenateSessions();
mysort.plot.SliderDataAxes(MSM)