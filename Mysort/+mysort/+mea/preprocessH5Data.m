function preprocessH5Data(MEA, varargin)
    P.smad_thr = 4;
    P.smad_Tf = 100;
    P.smad_maxLen = 500000;
    P.sd_thr = 3.5;
    P.sd_minPeakDistance = 20;
    P.sd_maxLen = 1000000;
    P.sd_chunkSize = 100000;
    P.cut_Tf = 30;
    P.cut_left = 5;
    P.cut_upsample = 8;
    P.outFileName = [];
    P = mysort.util.parseInputs(P, varargin, 'error');
    
    assert(isa(MEA, 'mysort.mea.CMOSMEA'), 'MEA must be a mysort.mea.CMOSMEA object');
    [Len nC] = size(MEA);
    if isempty(P.outFileName)
        P.outFileName = [MEA.fname '.preproc.mat'];
    end
    
    if exist(P.outFileName, 'file')
        warning('File already exists! Delete if necessary');
        return
    end
    
    % Calculate channel wise noise standard deviation with the median
    % absolute deviation (MAD), invert data to ignore negative peaks
    % for that calculation
    disp('Computing noise std...'); tic
    smadL = min(Len, P.smad_maxLen);
    smad = mysort.noise.estimateSigma(...
        -MEA(1:smadL,1:nC), P.smad_Tf, P.smad_thr);
    disp('Done.'); toc

    % Detect spikes in the beginning of the file
    disp('Detecting spikes...'); tic
%     sdL = min(Len, P.sd_maxLen);
    chunker = mysort.util.Chunker(Len, 'chunkSize', P.sd_chunkSize, 'progressDisplay', 'console', 'minChunkSize', 1000);
    mi_pks = {};  times = {};
%     ma_pks = {};
    for c=1:nC
        mi_pks{c,1} = [];
%         ma_pks{c,1} = [];
        times{c,1} = [];
    end
    while chunker.hasNextChunk()
        chunk = chunker.getNextChunk();
        X = double(MEA(chunk(1):chunk(2),1:nC));
        for c=1:nC
            [pks_, times_] = findpeaks(-X(:,c), 'MINPEAKHEIGHT', smad(c)*P.sd_thr,...
                                  'MINPEAKDISTANCE', P.sd_minPeakDistance);
            pks_ = pks_(:);
            times_ = times_(:);
            mi_pks{c,1} = [mi_pks{c}; pks_];
            times{c,1} = [times{c}; times_+chunk(1)-1];
        end
    end
    disp('Done.'); toc
    
%     % Cut spikes and upsample
%     disp('Cutting and upsampling spikes...'); tic
%     unlimited = H5ML.get_constant_value('H5S_UNLIMITED');
%     dims = [1 P.cut_Tf*P.cut_upsample];
%     maxDims = [unlimited P.cut_Tf*P.cut_upsample];
%     h5type = 'H5T_NATIVE_DOUBLE';
%     chunk_dims = dims;
%     deflation = 0;
%     spike_fname = [P.outFileName '_spks.h5'];
%     spike_h5path = '/cutspikes';
%     SP = mysort.h5.matrixCreate(spike_fname, spike_h5path, dims, maxDims,...
%         h5type, chunk_dims, deflation);
%     chan_time_amps_h5path = '/chan_time_amps';
%     CHAN_TIME_AMPS = mysort.h5.matrixCreate(spike_fname, chan_time_amps_h5path,...
%         [1 5], [unlimited 5], h5type, [1000 5], deflation);
%     sidx = 1;
%     for c=1:nC
%         for s=1:length(mi_pks{c,1})
%             s1 = max(mi_pks{c,1}(s)-P.cut_left, 1);
%             s2 = min(mi_pks{c,1}(s)-P.cut_left+P.cut_Tf-1, size(MEA,1));
%             sp = double(MEA(s1:s2,c));
%             rsp = resample(sp, P.cut_upsample, 1);
%             [mi mi_idx] = min(rsp);
%             [ma ma_idx] = max(rsp);
%             SP(sidx, 1:length(rsp)) = rsp';
%             real_mi_time = s1+(mi_idx-1)/P.cut_upsample;
%             real_ma_time = s1+(ma_idx-1)/P.cut_upsample;
%             mi_pks{c,1}(s) = mi;
%             ma_pks{c,1}(s) = ma;
%             times{c,1}(s) = real_mi_time;
%             CHAN_TIME_AMPS(sidx, :) = [c real_mi_time mi real_ma_time ma]; 
%             sidx=sidx+1;
%         end
%     end
%     disp('Done.'); toc

    % Build sorted matrix containing time points, channel index and
    % amplitude of found spikes
    disp('Building spike index...'); tic
    times(cellfun(@isempty, times)) = [];
    mi_pks(cellfun(@isempty, mi_pks)) = [];
    allTimes = cell2mat(times);
    allAmps   = cell2mat(mi_pks);
    allChans  = zeros(size(allTimes));
    nCum = 0;
    for i=1:length(times)
        allChans(nCum+1:nCum+length(times{i}),1) = i;
        nCum = nCum + length(times{i});
    end
    timesChansAmpsSorted = sortrows([allTimes allChans allAmps]);
    fprintf('maxtime: %d\n', max(allTimes));
    disp('Done.'); toc
    
    % Save the data
    version = 1;
    readme = 'This file was created by mysort.mea.prepareH5.m and contains preprocessed information about the associated h5 mea recording. Do not change if you want the mea objects to work properly with it';
    date_ = date();
    save(P.outFileName, 'timesChansAmpsSorted', 'P', 'smad', 'version', ...
        'readme', 'date_', '-v7.3');