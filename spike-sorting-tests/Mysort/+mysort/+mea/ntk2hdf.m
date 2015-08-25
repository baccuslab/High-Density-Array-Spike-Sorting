function ntk2hdf(fpath, ntkfiles, hdfile, varargin)
    P.deflation = 1;
    P.type = 'float';
    P.chunkSize = 200;
    P.writelog = false;
    P = mysort.util.parseInputs(P, varargin, 'error');
    
    if ~iscell(ntkfiles)
        ntkfiles = {ntkfiles};
    end
    
    % check if this is linux
    if isempty(strfind(computer, 'LNX'))
        fprintf('This function only works under linux with the ntk2hdf converter installed!');
        return
    end        
    
    % check if outfile already exists
    if exist(hdfile, 'file')
        fprintf('Warning, output h5 file (%s) already exists. Appending currently not possible, remove or change filename.\n', hdfile);
        return
    end
    
    % check if all files exist:
    for i=1:length(ntkfiles)
        ffile = fullfile(fpath, ntkfiles{i});
        if ~exist(ffile, 'file')
            fprintf('Warning, could not find file: %s\n\nOperation terminated, no file created.\n', ffile);
            return
        end
    end
% [frankef@bs-dsvr14 ~]$ ntk2hdf -?
% 
% ************ Convert NTK to HDF5 format ************
% 
% ntk2hdf: unrecognized option `-?'
% the message as passed by cmdline is: -
% usage: ntk2hdf [--chunkSize frames] [--defaltion 0-9] [--output outputName] tracefile.ntk [ additional.trace.files ]
% 
%  * help       (h)  print this help
%  * chunkSize  (c)  Internally HDF5 organizes the file in this chunkSize
%  * deflation  (d)  use deflation 0-9 (default: 6)
%  * output     (o)  Use this filename for the converted .h5 file
%  * version         A directory which should be analyzed file by file
% 
% 
% tool version: 1.1
    logstr = '';
    if P.writelog
        [h5path h5name h5extension] = fileparts(hdfile);
        logFile = [h5path h5name '.ntkconversion.log'];
        logstr = [' 2> ' logFile];
        fprintf('Logging to %s\n', logFile);
    end
    % create call string for 
    str = ['! ntk2hdf --chunkSize ' num2str(P.chunkSize) ' -d ' num2str(P.deflation) ' -o ' hdfile];
    for i=1:length(ntkfiles)
        ffile = fullfile(fpath, ntkfiles{i});
        str = [str ' ' ffile ];
    end
    str = [str logstr];
    fprintf('starting conversion in terminal:\n %s\n', str);
    eval(str);
    fprintf('Done.\n');
    
    
    
    