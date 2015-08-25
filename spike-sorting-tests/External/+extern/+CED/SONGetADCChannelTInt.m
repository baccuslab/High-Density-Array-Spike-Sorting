function[data,flag,h,Tend]=SONGetADCChannelTInt(fid,chan,varargin)
% Reads an ADC (waveform) channel from a SON file.
%
%

% Malcolm Lidierth 02/02
% modified by Clemens Boucsein 02/05: now with the option to specify the
% time interval for which ADC data are extracted (varargin{1} = start time,
% varargin{2} = end time in seconds)
% 
% 060312: the last read block is returned, and can be fed back to the function
% at the next call as varargin(3). If present, read out starts at the block 
% next to the last read block, irrespective of the specified starting time
import extern.CED.*;
flag = 0;

FileH=SONFileHeader(fid);
SizeOfHeader=20;                                            % Block header is 20 bytes long
header=SONGetBlockHeaders(fid,chan);
Info=SONChannelInfo(fid,chan);
maxFileTime=SONTicksToSeconds(fid,Info.maxChanTime);
if(Info.kind==0) 
   warning('SONGetADCChannelTInt: No data on that channel');
   return;
end;

%**************************Sample Interval************************
switch Info.kind                                                    % Disk block headers
case {1,6,7,9}
    switch FileH.systemID
    case {1,2,3,4,5}                                                % Before version 6
        if (isfield(Info,'divide'))
           Interval=Info.divide*FileH.timePerADC;                  % Convert to microseconds
           TicsInterval = Info.divide;
           h.divide = Info.divide;
        end;
        
    case {6}                                                        % Version 6
        Interval=Info.lChanDvd*FileH.timePerADC;
        TicsInterval = Info.lChanDvd;
        h.divide = Info.lChanDvd;
    end;
otherwise
    warning('SONGetSampleInterval: Invalid channel type');
    return
end;



if length(varargin) > 1
   startTime = varargin{1};
   endTime = varargin{2};
else
   startTime = [];
   endTime = [];
end

%*****************Convert Seconds to Clock Ticks**************************

if (isempty(startTime)) || (startTime == 0) 
  Tstart = 0;
else 
  Tstart = double((startTime/FileH.usPerTime)/FileH.dTimeBase);   % start time in clock tics
%   Tstart = double((startTime/(FileH.usPerTime*FileH.timePerADC))/FileH.dTimeBase);   % start time in clock tics
%   Tstart = double(startTime / (FileH.dTimeBase));   % start time in clock tics
%   Tstart = double(startTime / (FileH.usPerTime*FileH.timePerADC*FileH.dTimeBase));   % start time in clock tics
end

if length(varargin) > 2
    prevTend = varargin{3};
    Tstart = prevTend + TicsInterval;
%     Tstart = prevTend;
end

if startTime > maxFileTime
    disp('start time is bigger than file length... aborting');
    return
end

if (endTime > maxFileTime) || (isempty(endTime))
    disp('end time is too big! Using maximum of available data...');
   Tend = double((maxFileTime/FileH.usPerTime)/FileH.dTimeBase);
%    Tend = double((maxFileTime/(FileH.usPerTime*FileH.timePerADC))/FileH.dTimeBase);
%    Tend = double(maxFileTime/(FileH.dTimeBase));
%    Tend = double(maxFileTime/(FileH.usPerTime*FileH.timePerADC*FileH.dTimeBase));
   flag = 1;
else
   Tend = double((endTime/FileH.usPerTime)/FileH.dTimeBase);
%    Tend = double((endTime/(FileH.usPerTime*FileH.timePerADC))/FileH.dTimeBase);      % end time in clock tics
%    Tend = double(endTime / (FileH.dTimeBase));      % end time in clock tics
%    Tend = double(endTime / (FileH.usPerTime*FileH.timePerADC*FileH.dTimeBase));      % end time in clock tics
end

%************Find which block Tstart and Tend are in*************

startBlock = find (header(3,:)<= Tstart);
endBlock = find (header(3,:)<= Tend);


if isempty(startBlock)
  startBlock = 1;
else startBlock = max (startBlock) + 1;
end

if isempty(endBlock)
  endBlock = 1;
else endBlock = max (endBlock) + 1;
end

if flag == 1
    endBlock = size(header,2);
end


%*****************************************************************

SampleInterval=(header(3,1)-header(2,1))/(header(5,1)-1);   % Sample interval in clock ticks
NumberOfSamples = sum(header(5,startBlock:endBlock));
h.start=header(2,2);                                        % Time of first data point (clock ticks)
h.stop=header(3,2);                                         % End of data (clock ticks)

% keyboard;


NumFrames=1;                                                % Number of frames. Initialize to one.
Frame(1)=1;
for m=1:Info.blocks-1                                       % Check for discontinuities in data record
   IntervalBetweenBlocks=header(2,m+1)-header(3,m);
   if IntervalBetweenBlocks>SampleInterval                 % If true data is discontinuous (triggered)
       NumFrames=NumFrames+1;                              % Count discontinuities (NumFrames)
       Frame(m+1)=NumFrames;                               % Record the frame number that each block belongs to
   else
       Frame(m+1)=Frame(m);                                % Pad between discontinuities
   end;
end;

if NumFrames==1                                            % Continuous sampling - one frame only
   data=int16(zeros(1,NumberOfSamples));                   % Pre-allocate memory for data
   pointer=1;
   offsetWords = 0;
   for m = startBlock:endBlock
          offsetBytes = 2 * ceil((Tstart - header(2,m))/(Interval));
          if offsetBytes < 0                               % returns reading offset in bytes for the
             offsetBytes = 0;                              % first data block
          else
             offsetWords = offsetBytes/2;                  % store offset for calculation of proper data vector size
          end
          cutEndWords = floor((header(3,m) - Tend)/(Interval));
          if cutEndWords < 0                               % returns length of data stretch which should not
             cutEndWords = 0;                              % be read from the last block
          end
          lengthWords = header(5,m) - offsetBytes/2 - cutEndWords;
          fseek(fid, header(1,m) + offsetBytes + SizeOfHeader,'bof');
          data(pointer:pointer + lengthWords-1) = fread(fid, lengthWords ,'int16=>int16');
          pointer = pointer + lengthWords;
   end;
   dataLength = sum(header(5,startBlock:endBlock)) - offsetWords - cutEndWords;
   data = data(1:dataLength);   
else                                                       % Frame based data -  multiple frames
   FrameLength=NumberOfSamples/NumFrames;                  % Data points to a frame - assumed constant for all frames
   data=int16(zeros(NumFrames,FrameLength));               % Pre-allocate array
   start=1;                                                % Pointer into array for each disk data block
   Frame(Info.blocks+1)=-99;                               % Dummy entry to avoid index error in for loop
   for m=1:Info.blocks                                        
       fseek(fid,header(1,m)+SizeOfHeader,'bof');
       data(Frame(m),start:start+header(5,m)-1)=fread(fid,header(5,m),'int16=>int16');
       if Frame(m+1)==Frame(m)
           start=start+header(5,m);                        % Increment pointer or.....
       else
           start=1;                                        % begin new frame
           h.start(Frame(m))=header(2,m);                  % Time of first data point in frame (clock ticks)
           h.stop(Frame(m))=header(3,m);                   % End time, clock ticks
       end;
   end;
end;


h.start=SONTicksToSeconds(fid,h.start);
h.stop=SONTicksToSeconds(fid,h.stop);

if(nargout>1)
h.FileName=Info.FileName;                                   % Set up the header information to return
h.system=['SON' num2str(FileH.systemID)];                   % if wanted
h.FileChannel=chan;
h.phyChan=Info.phyChan;
h.kind=Info.kind;
h.timePerADC = FileH.timePerADC;
% h.lChanDvd = Info.lChanDvd;
% h.divide = Info.divide;
h.usPerTime = FileH.usPerTime;
h.dTimeBase = FileH.dTimeBase;

h.blocks=Info.blocks;
h.preTrig=Info.preTrig;

h.comment=Info.comment;
h.title=Info.title;
h.sampleinterval=SONGetSampleInterval(fid,chan);
h.scale=Info.scale;
h.offset=Info.offset;
h.units=Info.units;
end;

 
  
  
             
