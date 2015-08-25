function prefilterH5Data(ffile, varargin)
    P.hpf = 400;
    P.lpf = 6000;
    P.filterOrder = 10;
    P.filter_chunkSize = 50000;
    P.filtfilt = 0;
    P.deflation = 1;
    P.h5chunkLen = 1000;
    P.outFileName = '';
    P.reverseX = 0;
    P = mysort.util.parseInputs(P, varargin, 'error');

    if isempty(P.outFileName)
        P.outFileName = [ffile(1:end-2) 'filtered' ...
                         num2str(P.filterOrder) '_' num2str(P.hpf) ...
                            '_' num2str(P.lpf) '.h5'];
    end
    if exist(P.outFileName, 'file')
        warning('The prefiltered file already exists. Delete if necessary');
        return
    end
    assert(P.filterOrder/2 == round(P.filterOrder/2), 'filterOrder must be multiple of 2');
    MEA = mysort.mea.CMOSMEA(ffile);

    T = size(MEA,1);
    nC = size(MEA,2);

    HD = mysort.mea.filter_design(P.hpf, P.lpf, MEA.samplesPerSecond, P.filterOrder);
    HD.PersistentMemory = 1;
    % init filter to reduce filter artifact at beginning
    R = MEA.getScaledData(1:1000,:);   
    Y = filter(HD, R);
    
    h5path = '/prefiltered';
    dims = size(MEA);
    maxDims = size(MEA);
    h5type = 'H5T_NATIVE_FLOAT';
    chunkDims = [P.h5chunkLen nC];
    X = mysort.h5.matrixCreate(P.outFileName, h5path, dims, maxDims, h5type, chunkDims, P.deflation);
    
    fprintf('Starting forward filtering...\n');
    chunker = mysort.util.Chunker(T, 'chunkSize', P.filter_chunkSize,...
        'progressDisplay', 'console', 'chunkOverlap', 0, 'minChunkSize', 1000);
    while chunker.hasNextChunk()
        chunk = chunker.getNextChunk();
        R = MEA.getScaledData(chunk(1):chunk(2),:);
        Y = filter(HD, R);
        Y = single(Y);
        X(chunk(1):chunk(2),:) = Y;
    end
    
    
    readme = 'This file was created by mysort.mea.prefilterH5Data.m and contains preprocessed information about the associated h5 mea recording. Do not change if you want the mea objects to work properly with it';
    hdf5write(P.outFileName, '/P', P, '/date', date, '/version', 1, '/readme', readme, 'WriteMode','append');
    
    if ~P.filtfilt;
        % suppress filter artifact
        %X(1:100, :) = 0;
        return
    end
    
    h5path = '/prefiltfiltered';
    Z = mysort.h5.matrixCreate(P.outFileName, h5path, dims, maxDims, h5type, chunkDims, P.deflation);

    % init filter to reduce filter artifact at beginning
    R = flipud(X(end-1000:end,:));   
    Y = filter(HD, R);
    
    fprintf('Starting backward filtering...\n');
    chunker = mysort.util.Chunker(T, 'chunkSize', P.filter_chunkSize,...
        'progressDisplay', 'console', 'chunkOverlap', 0, 'minChunkSize', 1000);
    while chunker.hasNextChunk()
        chunk = chunker.getNextChunk();
        a=T-chunk(2)+1;
        b=T-chunk(1)+1;
        R = X(a:b,:);
        R = flipud(R);
        Y = filter(HD, R);
        Y = flipud(Y);
        Y = single(Y);
        Z(a:b,:) = Y;
    end
    
    % suppress filter artifact
%     Z(T-150:T, :) = 0;
%     Z(1:150, :)   = 0;
%     X(1:100, :)   = 0;
%     