spsort_path = '/links/groups/hima/recordings/HiDens/SpikeSorting/';
dpath = '/net/bs-filesvr01/export/group/hierlemann/recordings/HiDens/'; %~/bel.svn/hima_internal/cmosmea_recordings/trunk/
configfile_path = '/home/frankef/bel.svn/hima_internal/cmosmea_recordings/trunk/';


enames = {'19Apr2012_DSGCs_rabbit'
    '19Jun2012_DSGCs_rabbit'
    '23Aug2012_DSGCs_rabbit'
    '30Aug2012_DSGCs_rabbit'
    '18Sep2012_DSGCs_rabbit'
    '20Sep2012_DSGCs_rabbit'
    '8Nov2012_DSGCs_rabbit'
    '27Nov2012_DSGCs_rabbit'
    '15Jan2013_DSGCs_rabbit'
    '17Jan2013_DSGCs_rabbit'
    '22Jan2013_DSGCs_rabbit'
    '24Jan2013_DSGCs_rabbit'};

egroup = 'Roska';
spikeSortingRunName = 'run_felix5';
nCPUs = 1;

for expi = 3%1:1%length(enames)
    % set name of experiment
    ename = enames{expi};
    % load configuration file
    cfile = fullfile(configfile_path, egroup, ename, 'Matlab', 'configurations.mat');
    C = load(cfile);
    % this is the active config that will be sorted
    config_idx = 1;
    
    ntkfilelist = C.flist(C.configs(1).flistidx);
    % set some paths
    egroup_dpath = fullfile(dpath, egroup);
    exp_dpath = fullfile(egroup_dpath, ename);
    ntk_dpath = fullfile(exp_dpath, 'proc');
    spsort_out_path = fullfile(spsort_path, egroup, ename);

    % This container will store the template chains. For each chain it says
    % which template ID in which ntk file belongs to this chain.
    R = struct();
    TChains = {};
    % convert to h5
    for i=1:1%length(ntkfilelist)
        % Convert ntk 2 hdf5
        [PATH,NAME,EXT] = fileparts(ntkfilelist{i});
        ntkfile = [NAME EXT];
        hd5file = fullfile(spsort_out_path, [NAME '.h5']);
        gdffile = fullfile(spsort_out_path, [NAME '.gdf']);
        spath = fullfile(spsort_out_path, NAME, 'sortings', spikeSortingRunName);
        groupfile = fullfile(spath, [spikeSortingRunName '_el_groups.mat']);
        
        load(groupfile);
        for g = 1:length(groups)
            umsfile = fullfile(spath, [spikeSortingRunName sprintf('_group%03d', g) '_Export4UMS2000.mat']);
            UMS = load(umsfile);
            
        end
        
        sorting_result_file = fullfile(spath, sprintf('group%03d', i), 'run_felix5_templates.mat');
        
        T = load(sorting_result_file);
        if i==1
            % In the first run, open only new TChains
        else
            % find pairs
            

            
            % for those who are paired put them to the respective TChain
            
            % the others start a new TChain
        end        
    end
    
    final_result_file = fullfile(spsort_out_path, [spikeSortingRunName '_sorting']);
end

% %% Plot electrode groupings
% figure;
% ah1= axes();
% for i=1:length(groups)
%     x = ME.electrodePositions(groupsidx{i},1);
%     y = ME.electrodePositions(groupsidx{i},2);
%     plot(ah1, x+1, y, 'x', 'color', mysort.plot.vectorColor(i), 'markersize', 14, 'linewidth', 2);
%     hold on
% end
% % linkaxes(ah, 'xy');
% figure;
% ah = mysort.plot.subplots(length(groups));
% for i=1:length(groups)
%     x = ME.electrodePositions(groupsidx{i},1);
%     y = ME.electrodePositions(groupsidx{i},2);
%     plot(ah(i), x+1, y, 'x', 'color', mysort.plot.vectorColor(i), 'markersize', 14, 'linewidth', 2);
%     hold on
% end
% linkaxes(ah, 'xy');
% set(ah, 'xlim', get(ah1, 'xlim'), 'ylim', get(ah1, 'ylim'));

