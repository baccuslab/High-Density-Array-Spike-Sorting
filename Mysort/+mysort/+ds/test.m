gME = self.wfManager.MultiElectrode;
gENr = gME.electrodeNumbers;

[T1, ME1, wfsAllSessions1] = self.wfManager.getTemplate4Idx(idx, cutLeft, cutLength, gENr(1:2));
[T2, ME2, wfsAllSessions2] = self.wfManager.getTemplate4Idx(idx, cutLeft, cutLength, gENr(1:4));
[T3, ME3, wfsAllSessions3] = self.wfManager.getTemplate4Idx(idx, cutLeft, cutLength, gENr(12));

figure; 
subplot(3,1,1)
plot(T1');

subplot(3,1,2)
plot(T2');

subplot(3,1,3)
plot(T');