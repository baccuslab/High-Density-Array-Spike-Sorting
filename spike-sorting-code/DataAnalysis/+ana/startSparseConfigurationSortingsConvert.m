spsort_path = '/links/groups/hima/recordings/HiDens/SpikeSorting/';
dpath = '/net/bs-filesvr01/export/group/hierlemann/recordings/HiDens/'; %~/bel.svn/hima_internal/cmosmea_recordings/trunk/
configfile_path = '/home/frankef/bel.svn/hima_internal/cmosmea_recordings/trunk/';


enames = {'19Apr2012_DSGCs_rabbit'
    '19Jun2012_DSGCs_rabbit'
    '23Aug2012_DSGCs_rabbit'
    '30Aug2012_DSGCs_rabbit'
    '18Sep2012_DSGCs_rabbit'
    '20Sep2012_DSGCs_rabbit' % spelling mistake in network folder! not svn!
    '8Nov2012_DSGCs_rabbit'
    '27Nov2012_DSGCs_rabbit'
    '15Jan2013_DSGCs_rabbit'
    '17Jan2013_DSGCs_rabbit'
    '22Jan2013_DSGCs_rabbit'
    '24Jan2013_DSGCs_rabbit'};

egroup = 'Roska';
for expi = 6:length(enames)
    % set name of experiment
    ename = enames{expi};
    % load configuration file
    cfile = fullfile(configfile_path, egroup, ename, 'Matlab', 'configurations.mat');
    C = load(cfile);
    % this is the active config that will be sorted
    config_idx = 1;
    spikeSortingRunName = 'run_felix1';
    ntkfilelist = C.flist(C.configs(1).flistidx);
    % set some paths
    group_dpath = fullfile(dpath, egroup);
    exp_dpath = fullfile(group_dpath, ename);
    % treat the spelling mistake in the 20Sep experiment!
    if ~isempty(strfind('20Sep2012_DSGCs_rabbit', ename))
        ntk_dpath = fullfile(fullfile(group_dpath, '20ep2012_DSGCs_rabbit'), 'proc');
    else
        ntk_dpath = fullfile(exp_dpath, 'proc');
    end
    spsort_out_path = fullfile(spsort_path, egroup, ename);
    
    % convert to h5
    for i=1:length(ntkfilelist)
        % Convert ntk 2 hdf5
        [PATH,NAME,EXT] = fileparts(ntkfilelist{i});
        ntkfile = [NAME EXT];
        hd5file = fullfile(spsort_out_path, [NAME '.h5']);
        gdffile = fullfile(spsort_out_path, [NAME '.gdf']);
        spath = fullfile(spsort_out_path, NAME, 'sortings');
        if ~exist(spath, 'file')
            mkdir(spath);
        end    
        mysort.mea.ntk2hdf(ntk_dpath, ntkfile, hd5file, 'writelog', true);

        % build mea object
        sessionmea = mysort.mea.CMOSMEA(hd5file, 'hpf', 350, 'lpf', 6000, 'filterOrder', 6);
        % prefilter
        sessionmea.prefilter();    
    end
end

