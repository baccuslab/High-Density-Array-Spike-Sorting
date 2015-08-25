function chanList = ListChannels(filename)
import extern.CED.*;
fileInfo = GetFileInformation(filename);
disp(['ChNr' '    ' 'title']);
for m = 1:length(fileInfo.ChanList)
    List{m} = [num2str(fileInfo.ChanList(m).number) '     ' fileInfo.ChanList(m).title];
    disp(List{m});
end

fid = fopen(filename);
Info=SONChannelInfo(fid,fileInfo.ChanList(1).number);
fclose(fid);
maxTime = Info.maxChanTime/100000;
disp(['max. Time: ' num2str(maxTime) ' seconds']);