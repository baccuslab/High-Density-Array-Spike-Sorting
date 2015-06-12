
function[FileH,Info,data]=testSONGetADCChanTInt()
import extern.CED.*;
[filename,pathname] = uigetfile ('*.smr','Select .smr data file');
fullFilename = strcat(pathname, filename);

startTime=10;
endTime=20;
VoltChan=2;
fid=fopen(fullFilename);

[data,flag] = SONGetADCChannelTInt(fid,VoltChan,startTime,endTime);
Info=SONChannelInfo(fid,VoltChan);
FileH=SONFileHeader(fid);

fclose(fid);

FileInfo = GetFileInformation(fullFilename);
SampleFreqVolt = 1/FileInfo.ChanList(VoltChan).sampleInt; 

dataLength=length(data)*FileInfo.ChanList(VoltChan).sampleInt;
disp(['filename: ' filename]);
disp(['cut out data length should be: ' num2str(endTime-startTime) ' sec']);
disp(['cut out data length actually is: ' num2str(dataLength) ' sec']);
disp(['us per time: ' num2str(FileH.usPerTime)]);
disp(['time per ADC: ' num2str(FileH.timePerADC)]);
