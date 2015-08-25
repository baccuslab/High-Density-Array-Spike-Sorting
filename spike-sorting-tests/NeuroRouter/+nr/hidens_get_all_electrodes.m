function els = hidens_get_all_electrodes(version)
    packagename = 'nr';
    bufferfname = 'hidens_get_all_electrodes.mat';
    s = what(packagename);
    
    datafile = fullfile(s.path, bufferfname);
    D= load(datafile);
    if version == 1
        els = D.els1;
    elseif version == 2
        els = D.els2;
    else
        error('version can be either 1 or 2!');
    end
        