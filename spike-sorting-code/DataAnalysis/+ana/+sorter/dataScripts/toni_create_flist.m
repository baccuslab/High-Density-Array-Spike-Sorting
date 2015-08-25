function configFile = toni_create_flist(expname, configPath)
    pd = pdefs();
    % make sure you have the path of a folder with all .ntk files WITH A SINGLE
    % COMMON configuration in the variable "confPath"
    % expname = 'Cortex mice';

    expDataPath          = fullfile(pd.mea1kData, 'Antonia', expname);
    
    configNTKLists = {};
    xy_pos = {};

    configFile = fullfile(configPath, [expname '.mat']);
    
    if exist(configFile, 'file')
        disp('Config File exists already! Loading !');
        return
    end

    % NEW VERSION FOR MEA 1K
    h5Files = dir(fullfile(expDataPath, '*.h5'));
    nrkFiles = dir(fullfile(expDataPath, 'config', '*.mapping.nrk'));
    assert(length(nrkFiles) == 1, 'Must be exactly one .nrk file!');

    xy_pos{1} = fullfile(expDataPath, 'config', nrkFiles(1).name);

    for i=1:length(h5Files)
        configNTKLists{1}{i} = fullfile(expDataPath, h5Files(i).name);
    end

    %%
    if ~exist(configPath, 'file'); mkdir(configPath); end
    save(configFile, 'configNTKLists', 'xy_pos');
end