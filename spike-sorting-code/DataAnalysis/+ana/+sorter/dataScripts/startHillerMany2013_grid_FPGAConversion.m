% matlabpool(12);

expNames = {'configs_19Mar2015_DZ_mea1k'}; 

%% PREPROCESSING AND SUBMISSION TO GRID ENGINE
all_Experiments_H5Files = {};
for eidx = 1:length(expNames)
    expName = expNames{eidx};
    fprintf('Preprocessing %d\n', expName);
    
    %% INIT
    pd = pdefs();
    inputPath = fullfile(pd.networkTempShare, 'Hillier_2013');
    project_name = [expName '_proj'];
    outSuffix = 'Out';


    %% PREPROCESS MEA1K FILES
    outPath = fullfile(inputPath, [expName outSuffix]);
    if ~exist(outPath, 'file'); mkdir(outPath); end
    confListFile = fullfile(inputPath, [expName '.mat']);
    c = load(confListFile);
    configNTKLists = c.configNTKLists;
    bIsMea1Krecording = true;
    try
        nrk_files_per_config = c.xy_pos;
    catch
        % THIS IS AN OLD HIDENS RECORDING
        bIsMea1Krecording = false;
        nrk_files_per_config = [];
    end

    %% Now we CONVERT all ntk files to H5. Doing this we also prefilter the data
    allFiles = {};
    allCombis = [];
    h5FileLists = {};
    mapPerFile = {};
    for i=1:length(configNTKLists)
        if bIsMea1Krecording
            % load config file
            fname = nrk_files_per_config{i};
            if iscell(fname)
                fname = fname{1};
            end            
            fname = strrep(fname, '.u.', '.');
            try 
                X = dlmread(fname);
            catch
                fname = strrep(fname, 'config/config_', 'config/configs_');
                try
                    X = dlmread(fname);
                catch
                    fname = strrep(fname, 'config/config', 'config/config_');
                    X = dlmread(fname);
                end
            end
            isNotConnected = any(X==-1,2);
            X(isNotConnected,:) = [];
            % Im H5 file ist spalte 0 channel 0.
            % 1023 ist der letzte Kanal dann DAC und !command counter"
            % DAC ist nicht im nrk file 
            idx = X(:,1) + 1;
            m = struct();
            m.mposx(idx) = X(:,2);
            m.mposy(idx) = X(:,3);
            m.elNo(idx)   = X(:,1);
            m.chans(idx)  = X(:,1);
            m.noChans = length(m.elNo);
        end
        for k=1:length(configNTKLists{i})
            fn = configNTKLists{i}{k};
            allFiles{end+1} = fn;
            if bIsMea1Krecording
                mapPerFile{end+1} = m;
            end
            [ffolder ffname fsuffix] = fileparts(fn);
            fnpf = fullfile(outPath, [ffname fsuffix '_prefilterFPGA.h5']);
            h5FileLists{i}{k} = fnpf;
        end
        all_Experiments_H5Files{eidx} = h5FileLists;
    end

    tic
    for combi=1:length(allFiles)
        disp(combi)
        fn = allFiles{combi};
        [ffolder ffname fsuffix] = fileparts(fn);
        fnpf = fullfile(outPath, [ffname fsuffix '_prefilterFPGA.h5']);
        bConvertFile = true;
        disp(['Checking if file exists ' num2str(combi)])
        if exist(fnpf, 'file')
            try
                bIsInProcess = hdf5read(fnpf, '/bFileIsInProcess');
    %                 bIsInProcess = 1;
                if bIsInProcess
                    % This means the file is still in process. We assume no
                    % other user/matlab instance is writing so it is
                    % probably a left over from some other run and we
                    % delete it.
                    delete(fnpf); 
                else
                    % The file was already properly processed, we dont need
                    % to do anything.
                    bConvertFile = false;
                end
            catch
                delete(fnpf);
            end
        end
        if bConvertFile
            % Only convert file if it was non existent or deleted
            if bIsMea1Krecording
                disp(['Converting Mea1K ' num2str(combi)])
                mysort.mea.preprocessMea1kH5File(fn, mapPerFile{combi}, 'prefilter', 1, 'outFile', fnpf, 'bUseFPGA_IIR_Filter', true);
            else
                error('This script is not done for Hidens data!');
                disp(['Converting Highdens ' num2str(combi)])
                mysort.mea.convertNTK2HDF(fn, 'prefilter', 1, 'outFile', fnpf);
            end
        end
    end
    fn = [];
    toc
    disp('Done preprocessing!')
    
end
 