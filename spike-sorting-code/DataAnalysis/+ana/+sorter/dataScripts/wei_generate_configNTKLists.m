%cd /home/wgong/bel.svn/hima_internal/cmosmea_recordings/trunk/organotypic/2013_12_11_1281/proc; 

function generate_configNTKLists(special_dir) 
%  Updated from the GENERATE_SEPARATE_FLISTS 
%  GENERATE_SEPARATE_FLISTS  generate flists for recordings, separated for
%  individual recording sessions.
%
%  GENERATE_SEPARATE_FLISTS('activity_scan')  generate flist for files
%  recorded in '../proc/activity_scan'.
%
%  See also generate_flist

%% gen flist
    special_dir=pwd

%%
files = dir([special_dir 'Trace*.ntk']);
c = {files.name};
[~, order] = sort([files(:).datenum]);
c = c(order);
configNTKLists = cellfun(@(x) {fullfile(special_dir, x)}, c, 'uniformOutput', false);
save('2013_12_11_1281_configs.mat', 'configNTKLists'); % should be saved in folder: /home/wgong/Desktop/organotypic/2013_12_11_1281/Configs

  