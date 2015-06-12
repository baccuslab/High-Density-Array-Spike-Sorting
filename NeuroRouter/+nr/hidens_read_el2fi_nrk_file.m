function [elidx amps] = hidens_read_el2fi_nrk_file(configPath, fname)
    ffname = fullfile(configPath, [fname '.el2fi.nrk2']);
    fid=fopen(ffname);
    assert(fid>0, 'Could not open el2fi.nrk2 file!');
    elidx=[];
    amps =[];
    tline = fgetl(fid);
    while ischar(tline)
        [tokens] = regexp(tline, 'el\((\d+)\).*\(([a-z][a-z]\d+), filter\)', 'tokens');%,filter
        elidx(end+1)=str2double(tokens{1}{1});
        amps(end+1) =routingString2Amplifier(tokens{1}{2});
        tline = fgetl(fid);
    end
    fclose(fid);
    
    function a = routingString2Amplifier(s)
        if strcmp(s(1:2), 'fr')
            offset = 0;
            a = offset + str2double(s(3:end));
        elseif strcmp(s(1:2), 'fb')
            offset = 62;
            a = offset - str2double(s(3:end));
        elseif strcmp(s(1:2), 'fl')
            offset = 98;
            a = offset - str2double(s(3:end));
        elseif strcmp(s(1:2), 'ft')
            offset = 99;
            a = offset + str2double(s(3:end));
        else
            error('unknown routing type');
        end
    end
end
