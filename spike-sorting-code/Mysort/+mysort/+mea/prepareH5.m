function prepareH5(fpath, fname, varargin)
    error('this function is obsolete, use the CMOSMEA object!');
    P.hpf = 400;
    P.lpf = 6000;
    P.filter_order = 10;
    P.smad_thr = 4;
    P.smad_Tf = 100;
    P.smad_maxLen = 500000;
    P.sd_thr = 4.5;
    P.sd_minPeakDistance = 20;
    P.sd_maxLen = 5000000;
    P.sd_chunkSize = 100000;
    P = mysort.util.parseInputs(P, varargin, 'error');

    ffile_in  = fullfile(fpath, fname);
    ffile_out = fullfile(fpath, [fname '_preproc.mat']);
    
    if ~exist(ffile_in, 'file')
        warning('Could not find specified H5 file!');
        return
    end
    if exist(ffile_out, 'file')
        warning('Preprocessed file already exist. Delete if necessary!');
        return
    end

    % Create filter buffer for H5 file
    disp('Prefiltering data...'); tic
    mysort.mea.prefilterH5Data(ffile_in, P);
    disp('Done.'); toc
    
    mysort.mea.preprocessH5Data();
    