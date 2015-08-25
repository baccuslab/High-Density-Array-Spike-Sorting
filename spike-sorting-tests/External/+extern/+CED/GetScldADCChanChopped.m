function [fileInfo,varargout]= GetScldADCChanChopped(filename,startTime,endTime,varargin)
import extern.CED.*;
% provides access to a variable number of channels of spike2 data files 
% (extension .smr). 
% 
% Expects the name (including path) of the data file (filename), 
% the start and end time of the data stretch that should be read out in seconds
% (startTime, endTime), and a variable list of channels from which data should 
% be read (varargin). 
% 
% Returns a structure (fileInfo) which contains information about the file, the
% number of channels and scaling factors ect.; varargout contains a variable 
% number of vectors each containing the start-end time stretch of data of one
% channel specified in 'varargin'.
% 
% Check for the available number of channels in a given data file and the channel
% number with seperate script 'ListChannels.m'
% 
% 

Chunksize = 0.5; % size of data chunks in seconds


% filename = 'D:\Dymphie\InVivo\2004\040603\040603_030000.SMR';
% filename = 'D:\temp\Sonja_Päckle_0501\040616_010000_350s.smr';
% filename = 'D:\Matlab\Test\040811_cortex1_2cell006_test.smr';
% filename = '\\Zero\CBexpData\LFP_Tetrode\2005\050426\050426_cortex1_1cell002.smr';
% filename = 'D:\Matlab\Test\041018_cortex1_1cell005_test.smr';

% [filename,pathname] = uigetfile ('*.smr','Select .smr data file');
% filename = strcat(pathname, filename);

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

resolution = 6553.6;        % resolution of the AD converter board, devided by 10 (here: 16 bit = 65536 steps)

loops = ceil((endTime - startTime)/Chunksize);

% keyboard;

for m = 1 : length(selectChan)
    chunkEnd = startTime + Chunksize;
    [rawData, flag, h, Tend] = SONGetADCChannelTInt(fid,selectChan{m},startTime,chunkEnd);    
    for k = 2 : loops
        chunkEnd = chunkEnd + Chunksize;
        if chunkEnd > endTime
            chunkEnd = endTime;
        end
        [newDataChunk, flag, h, Tend] = SONGetADCChannelTInt(fid,selectChan{m},startTime,chunkEnd,Tend);
        temp(k)=newDataChunk(1)==rawData(end);
        rawData = [rawData newDataChunk];
    end
    ChanInd = find([fileInfo.ChanList.number] == selectChan{m});
    varargout{m} = ([fileInfo.ChanList(ChanInd).scale] * (double(rawData(:))/resolution)) + fileInfo.ChanList(ChanInd).offset;
end

