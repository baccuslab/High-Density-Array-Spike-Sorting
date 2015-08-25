token_path = '/net/bs-gridfs/gridexport/scratch/frankef/tokens/';
log_file = '~/submit_demon.log';
cd(token_path)
pd = pdefs();
while 1
    fnames = dir(fullfile(token_path, 'start_*.tok'));
    if ~isempty(fnames)
        disp('Found token, processing...')
        for i=1:length(fnames)
            token_file = fullfile(token_path, fnames(i).name);
            T = load(token_file, '-mat');
            proj_name = T.jobName;
            exp_name  = T.expName;
            submit_token_file = fullfile(token_path, 'submitted', fnames(i).name);
            disp(token_file)
            disp(proj_name)
            
            sortPath = T.sortOutPath;

            submit_str = sprintf('sh ~/autosubmit.sh %s %s &', sortPath, proj_name);
            move_str   = sprintf('mv %s %s &', token_file, submit_token_file);
            %submit_str = sprintf('ls %s', proj_name);

            mysort.util.logToFile(log_file, submit_str)
            [status, result] = system(submit_str);
            mysort.util.logToFile(log_file, result)

            mysort.util.logToFile(log_file, move_str)
            [status, result] = system(move_str);
            mysort.util.logToFile(log_file, result)
            mysort.util.logToFile(log_file, '################################')
            pause(2.5);
            disp('Done processing. Waiting...')
%             pause(3600);
        end
    else
        % disp('Nothing found. Waiting...')
    end
%     cd(token_path)
%     return
    pause(5)
end

disp('Done.')