function [varargout]= TimeInterval(startTime,endTime,varargin)
import extern.CED.*;

% filename = 'D:\Dymphie\InVivo\2004\040603\040603_030000.SMR';
% filename = 'D:\temp\Sonja_Päckle_0501\040616_010000_350s.smr';
filename = 'D:\Matlab\Test\040811_cortex1_2cell006_test.smr';

% [filename,pathname] = uigetfile ('*.smr','Select .smr data file');
% filename = strcat(pathname, filename);

ADC = '12 bit';
ADC = '16 bit';

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
%         case '12 bit'
%             varargout{m} = ([fileInfo.ChanList(ChanInd).scale] * (double(rawData(:))-2049))/409.6;
        case '12 bit'
            varargout{m} = [fileInfo.ChanList(ChanInd).scale] * (((double(rawData(:))-2049))/409.6);
        otherwise
            disp('ADC bit depth not specified; no output variables will be generated!');
    end
end

