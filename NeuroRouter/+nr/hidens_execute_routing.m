function hidens_execute_routing(configPath, fname, specsFolderName)
    if nargin < 3
        specsFolderName = 'matlab_specs';
    end
    nr_exe='NeuroDishRouter';
    nrk_file = fullfile(configPath, specsFolderName, [fname '.neuropos.nrk']);
    out_file = fullfile(configPath, fname);
    
    str = sprintf('%s -n -v 2 -l %s -s %s\n', nr_exe, nrk_file, out_file);
    fprintf('running: %s', str);
    unix(str);
