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
    '24Jan2013_DSGCs_rabbit'
    '29Jan2013_DSGCs_rabbit'}1;

egroup = 'Roska';
spikeSortingRunName = 'run_felix5';
nCPUs = 0;

for expi = 4%length(enames)
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

    % convert to h5
    for i=1:length(ntkfilelist)
        % Convert ntk 2 hdf5
        [PATH,NAME,EXT] = fileparts(ntkfilelist{i});
        ntkfile = [NAME EXT];
        hd5file = fullfile(spsort_out_path, [NAME '.h5']);
        gdffile = fullfile(spsort_out_path, [NAME '.gdf']);
        spath = fullfile(spsort_out_path, NAME, 'sortings', spikeSortingRunName);
        groupfile = fullfile(spath, [spikeSortingRunName '_el_groups.mat']);
        if ~exist(spath, 'file')
            mkdir(spath);
        end    
        mysort.mea.ntk2hdf(ntk_dpath, ntkfile, hd5file, 'writelog', true);

        % build mea object
        sessionmea = mysort.mea.CMOSMEA(hd5file, 'hpf', 350, 'lpf', 6000, 'filterOrder', 6);
        % prefilter
        sessionmea.prefilter();    


        % build local electrode neighborhoods that will be sorted together
        ME = sessionmea.getMultiElectrode();
        if exist(groupfile, 'file')
            disp('###');
            disp('GROUP File already exists, using precomputed groups!');
            disp('###');
            load(groupfile);
        else
            [groupsidx nGroupsPerElectrode] = mysort.mea.constructLocalElectrodeGroups(ME.electrodePositions(:,1), ME.electrodePositions(:,2));
            % replace electrode indices with electrode numbers
            groups = {};
            for ii=1:length(groupsidx)
                groups{ii} = ME.electrodeNumbers(groupsidx{ii});
            end
            save(groupfile, 'groups', 'groupsidx', 'nGroupsPerElectrode');
        end
    

    %%     sort
        logfile = fullfile(spath, 'log.txt');
    %         try
            tic
            ana.sortConfig(spath, sessionmea, 1, groups, groupsidx, spikeSortingRunName, nCPUs, 0);
            toc
    %         catch
            mysort.util.logToFile(logfile, ['iteration: ' num2str(i)]);
            mysort.util.logLastErrToFile(logfile);
    %         end
    end
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

