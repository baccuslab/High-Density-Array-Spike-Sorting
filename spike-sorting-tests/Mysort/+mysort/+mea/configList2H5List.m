function h5lists = configList2H5List(path, expName, configList)
    if nargin < 3
        confListFile = fullfile(path, [expName '.mat']);
        c = load(confListFile);
        configList = c.configNTKLists;
    end
    h5lists = cell(size(configList));
    for i = 1:length(configList)
        for k=1:length(configList{i});
            [a, b, c] = fileparts(configList{i}{k});
            h5lists{i}{k} = fullfile(path, [expName 'Out'], [b c '_prefilter.h5']);
        end
    end