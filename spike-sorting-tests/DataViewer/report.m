
%%    INIT
DH = dataviewer.DataHandle('dbconfig', db.munk.DBconfig_Select, ...
    'log_function', @(x) fprintf(x));


path = 'C:/Data/Report/';


E = 'L011';
epath = fullfile(path, E);
if ~exist(epath, 'dir')
    mkdir(epath);
end

%% Meta Retrieval
blocks = DH.getBlocks('experiments', E);
tetrodes = DH.getTetrodes('experiments', E);

for b=1:size(blocks,1)
    bpath = fullfile(epath, blocks{b,2});
    if ~exist(epath, 'dir')
        mkdir(epath);
    end
    for tet = 1:size(tetrodes,1);
        tetpath = fullfile(bpath, sprintf('%02d',tetrodes{tet,2}));
        if ~exist(tetpath, 'dir')
            mkdir(tetpath);
        end
        
    end
end


