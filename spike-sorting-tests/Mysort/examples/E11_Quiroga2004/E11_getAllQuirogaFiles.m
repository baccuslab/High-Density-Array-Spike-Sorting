function files = E10_getAllQuirogaFiles(path_to_directory)

    fnames = dir([path_to_directory '*.mat']);
    files = {};
    count = 0;
    for i=1:length(fnames)
        if ~( ~any(strfind(fnames(i).name, 'short')) && ...
              ~any(strfind(fnames(i).name, 'Burst')) && ...
              ~any(strfind(fnames(i).name, 'times_')) && ...
              ~any(strfind(fnames(i).name, 'Test')) && ...
              ~any(strfind(fnames(i).name, 'Drift')) && ...
              ~fnames(i).isdir)
            continue
        end
        count = count +1;
        files{count} = fnames(i).name;
%         disp(fnames(i).name)
    end