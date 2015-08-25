function caller_name = getCallerName()
    caller_name = '';
    dbs = dbstack('-completenames');
    
    if ~isempty(dbs)
        [a b c] = fileparts(dbs(end).file);
        sepidx = strfind(a, filesep);
        if length(sepidx) < 1
            lastFolder = '';
        elseif length(sepidx) == 1 
            lastFolder = a(sepidx(end)+1:end);
        else
            lastFolder = a(sepidx(end-1)+1:end);
        end        
        if length(dbs) > 2
            caller_name = [dbs(end).name '() Line: '  num2str(dbs(end).line) ' in ' lastFolder ]; % dbs(end).name c
        else
            % 
            caller_name = 'unknown (shell?)';
        end
    end