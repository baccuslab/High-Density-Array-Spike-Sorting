pd = pdefs();

token_path = fullfile( pd.intermediateRoland, 'tokens' );
log_file = '~/submit_demon.log';
cd(token_path)

while 1
    fnames = dir(fullfile(token_path, 'start_*.mat'));
    if ~isempty(fnames)
        disp('Found token, processing...')
        for i=1:length(fnames)
            token_file = fullfile(token_path, fnames(i).name);
            t = load(token_file);
            submit_token_file = fullfile(token_path, 'submitted', fnames(i).name);
            
            disp(token_file)
            
            disp(['Sorting Location: ' t.sortingLocation])
            disp(['Sorting Name: ' t.sortingName])
            submit_str = sprintf('sh ~/autosubmit.sh %s %s', t.sortingLocation, t.sortingName);
            move_str   = sprintf('mv %s %s', token_file, submit_token_file);

            mysort.util.logToFile(log_file, submit_str)
            [status, result] = system(submit_str);
            mysort.util.logToFile(log_file, result)
            
            if status == 0
                disp('Submit successful')
                disp(result)
            else
                disp('Submit failed')
            end
            
            mysort.util.logToFile(log_file, move_str)
            [status, result] = system(move_str);
            mysort.util.logToFile(log_file, result)
            pause(2.5)
            disp('Done processing. Waiting...')
        end
    end
    pause(10)
end

disp('Done.')
