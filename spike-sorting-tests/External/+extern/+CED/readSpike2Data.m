function[Volt,Curr,Scal,Trig] = readSpike2Data(filename);
import extern.CED.*;
% reads in spike2 (*.smr) data and returns current,
% voltage and scaled signal data vectors 
% 
% data are stored in structures; the structure carries Channel number,
% sample interval and scaling factor. In addition, raw data vectors for
% voltage, current and scaled signal
% 
% uses the function GetFileName and SON-Toolbox

% C.Boucsein 050110

[fileInfo] = GetFileInformation(filename);

fid = fopen (filename);
if fid < 1  
    disp('error reading file');return        
end;

for i = 1 : length(fileInfo.ChanList)
    ChContent = fileInfo.ChanList(i).title;
    if strncmp(ChContent,'volt',4)
        Volt.ChanNr = fileInfo.ChanList(i).number;
        Volt.timestep = fileInfo.ChanList(i).sampleInt;
        Volt.scale = fileInfo.ChanList(i).scale;
        Volt.units = fileInfo.ChanList(i).units;
    elseif strncmp(ChContent,'curr',4)
        Curr.ChanNr = fileInfo.ChanList(i).number;
        Curr.timestep = fileInfo.ChanList(i).sampleInt;
        Curr.scale = fileInfo.ChanList(i).scale;
        Curr.units = fileInfo.ChanList(i).units;
    elseif strncmp(ChContent,'scld',4)
        Scal.ChanNr = fileInfo.ChanList(i).number;
        Scal.timestep = fileInfo.ChanList(i).sampleInt;
        Scal.scale = fileInfo.ChanList(i).scale;
        Scal.units = fileInfo.ChanList(i).units;
    elseif strncmp(ChContent,'Trig',4)
        Trig.ChanNr = fileInfo.ChanList(i).number;
        % Trig.timestep = fileInfo.ChanList(i).sampleInt;
    end
end

Volt.data = SONGetChannel (fid,Volt.ChanNr);
Curr.data = SONGetChannel (fid,Curr.ChanNr);

if exist('Scal')
    Scal.data = SONGetChannel (fid,Scal.ChanNr);
else
    Scal.data = [];
end

if exist('Trig')
    Trig.data = SONGetMarkerChannel (fid,Trig.ChanNr);
else
    Trig.data = [];
end


fclose(fid);

% Volt.sclData = (([Volt.scale] * double(Volt.data(:)))/6553.6);
% Curr.sclData = (([Curr.scale] * double(Curr.data(:)))/6553.6);
% Scal.sclData = (([Scal.scale] * double(Scal.data(:)))/6553.6);


