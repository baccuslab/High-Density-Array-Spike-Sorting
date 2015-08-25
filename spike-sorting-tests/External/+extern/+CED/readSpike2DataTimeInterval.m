function [varargout]= readSpike2DataTimeInterval(startTime,endTime,varargin)
import extern.CED.*;
filename = 'D:\projects\Variability\In_Vivo_Data_für_Martin\031023_cell2_001_sec166_311.smr';
% filename = 'D:\Dymphie\InVivo\2004\040603\040603_030000.SMR';
% filename = 'D:\temp\sliceI_gr1\5khz_2khzsamp.smr';
% filename = 'D:\temp\040920\040920_cortex1_1cell002.smr';

% [filename,pathname] = uigetfile ('*.smr','Select .smr data file');
% filename = strcat(pathname, filename);

ADC = '12 bit';                 % use in case of reading 12 bit MCRack data
% ADC = '16 bit';               % use in case of Spike2 data (16 bit encoded)

selectChan = varargin;

[fileInfo] = GetFileInformation(filename);

AvailChannels = [fileInfo.ChanList.number];

 for f = 1 : length(varargin)
     ChanCheck = find(varargin{f} == AvailChannels);
     if ~isempty(ChanCheck)
        break
     else disp(['Channel Number ' varargin{f}, 'does not exist']);return
     end
 end
 
        

fid = fopen(filename);
if fid < 1
    disp('error reading file'); return
end;

for m = 1 : length(selectChan)
    rawData = SONGetADCChannelTInt(fid,selectChan{m},startTime,endTime);
    ChanInd = find([fileInfo.ChanList.number] == selectChan{m});
    switch ADC
        case '16 bit'
            varargout{m} = ([fileInfo.ChanList(ChanInd).scale] * double(rawData(:)))/6553.6;
        case '12 bit'
            varargout{m} = ((double(rawData(:))/6553.6)*[fileInfo.ChanList(ChanInd).scale]) + fileInfo.ChanList(ChanInd).offset;
    end
end

