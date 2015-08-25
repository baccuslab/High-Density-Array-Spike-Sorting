function RES = matrixTestSpeedHelper(ffile, dims, h5type, chunkDims, deflation)
    ffilemat = [ffile '.mat'];
    h5path = 'foo';
    if exist(ffile, 'file'); delete(ffile); end
    if exist(ffilemat, 'file'); delete(ffilemat); end    
    R = randn(dims);
    pause(.1)
    maxDims = dims;
    
    RES = [];
    t = 1;
    
    tic
    M = mysort.h5.matrixCreate(ffile, h5path, dims, maxDims, h5type, chunkDims, deflation);   
    RES(t) = toc; t=t+1;

    tic
    M(:,:) = R; 
    RES(t) = toc; t=t+1;

    tic
    X = M(:,:);
    RES(t) = toc; t=t+1;

    tic 
    s1 = floor(rand*100)+1;
    s2 = floor(rand*10)+1;
    x = M(s1:s1+1000, s2:s2+10);
    RES(t) = toc; t=t+1;
    
    tic 
    save(ffilemat, 'X', '-v7.3');
    RES(t) = toc; t=t+1;

    clear X

    tic 
    load(ffilemat, 'X');
    RES(t) = toc; t=t+1;   
    delete(M)