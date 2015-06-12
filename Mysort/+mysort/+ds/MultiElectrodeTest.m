ePos = (-1:-10:-1000)';
eNum = 1:length(ePos);
sME = mysort.ds.MultiElectrode(ePos, eNum);
sME2 = sME.getSubElectrode(eNum(1:4));
sME2

%%
MSME = mysort.ds.MultiSessionMultiElectrode([sME sME2], [1 2]);
MSME

MSME2 = MSME.getSubElectrode(eNum(1:2))

MSME3 = MSME.getSubElectrode(eNum(3:6))

MSME4 = MSME.getSubElectrode(eNum(5:6))