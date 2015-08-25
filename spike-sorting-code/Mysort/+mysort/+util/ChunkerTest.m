% function ChunkerTest()
    %% TEST 1
    chunksize = 100;
    Len = 999;
    
    overlap = 20;
    
    chunksO =[  1 120
               81 220
              181 320
              281 420
              381 520
              481 620
              581 720
              681 820
              781 920
              881 999];
    chunks =[   1 100
              101 200
              201 300
              301 400
              401 500
              501 600
              601 700
              701 800
              801 900
              901 999];
    
    C = mysort.util.Chunker(Len, 'chunkSize', chunksize, ...
                                 'chunkOverlap', overlap);
	
	chunksO_ = zeros(size(chunks));
    chunks_  = chunksO;
    i = 1;
    while C.hasNextChunk()
        [chunksO_(i,:) chunks_(i,:)] = C.getNextChunk();
        i = i+1;
    end
    
    assert(~any(any(chunksO~=chunksO_)));
    assert(~any(any(chunks~=chunks_)));
    
    C.reset();
    
	chunksO_ = zeros(size(chunks));
    chunks_  = chunksO;
    i = 1;
    while C.hasNextChunk()
        [chunksO_(i,:) chunks_(i,:)] = C.getNextChunk();
        i = i+1;
    end
    
    assert(~any(any(chunksO~=chunksO_)));
    assert(~any(any(chunks~=chunks_)));
    
    
    
    

    %% TEST 2
    chunksize = 100;
    Len = 1001;
    minchunksize = 1;
    overlap = 20;
    
    chunksO =[  1 120
               81 220
              181 320
              281 420
              381 520
              481 620
              581 720
              681 820
              781 920
              881 1001
              981 1001];
    chunks =[   1 100
              101 200
              201 300
              301 400
              401 500
              501 600
              601 700
              701 800
              801 900
              901 1000
             1001 1001];
    
    C = mysort.util.Chunker(Len, 'chunkSize', chunksize, ...
                                 'chunkOverlap', overlap,...
                                 'minChunkSize', minchunksize);
	
	chunksO_ = zeros(size(chunks));
    chunks_  = chunksO;
    i = 1;
    while C.hasNextChunk()
        [chunksO_(i,:) chunks_(i,:)] = C.getNextChunk();
        i = i+1;
    end
    
    assert(~any(any(chunksO~=chunksO_)));
    assert(~any(any(chunks~=chunks_)));
    
    C.reset();
    
	chunksO_ = zeros(size(chunks));
    chunks_  = chunksO;
    i = 1;
    while C.hasNextChunk()
        [chunksO_(i,:) chunks_(i,:)] = C.getNextChunk();
        i = i+1;
    end
    
    assert(~any(any(chunksO~=chunksO_)));
    assert(~any(any(chunks~=chunks_)));
    
    
    %% TEST 3
    chunksize = 100;
    Len = 1001;
    minchunksize = 10;
    overlap = 20;
    
    chunksO =[  1 120
               81 220
              181 320
              281 420
              381 520
              481 620
              581 720
              681 820
              781 920
              881 1001
              972 1001];
    chunks =[   1 100
              101 200
              201 300
              301 400
              401 500
              501 600
              601 700
              701 800
              801 900
              901 991    % --> corrected len of 100-9
              992 1001]; % --> len of 10 !
    
    C = mysort.util.Chunker(Len, 'chunkSize', chunksize, ...
                                 'chunkOverlap', overlap,...
                                 'minChunkSize', minchunksize);
	
	chunksO_ = zeros(size(chunks));
    chunks_  = chunksO;
    i = 1;
    while C.hasNextChunk()
        [chunksO_(i,:) chunks_(i,:)] = C.getNextChunk();
        i = i+1;
    end
    
    assert(~any(any(chunksO~=chunksO_)));
    assert(~any(any(chunks~=chunks_)));
    
    C.reset();
    
	chunksO_ = zeros(size(chunks));
    chunks_  = chunksO;
    i = 1;
    while C.hasNextChunk()
        [chunksO_(i,:) chunks_(i,:)] = C.getNextChunk();
        i = i+1;
    end
    
    assert(~any(any(chunksO~=chunksO_)));
    assert(~any(any(chunks~=chunks_)));    
    
    
    %% TEST 4
    chunksize = 1000;
    Len = 1000;
    minchunksize = 10;
    overlap = 20;
    
    chunksO =[  1 1000];
    chunks =[   1 1000];
    
    C = mysort.util.Chunker(Len, 'chunkSize', chunksize, ...
                                 'chunkOverlap', overlap,...
                                 'minChunkSize', minchunksize);
	
	chunksO_ = zeros(size(chunks));
    chunks_  = chunksO;
    i = 1;
    while C.hasNextChunk()
        [chunksO_(i,:) chunks_(i,:)] = C.getNextChunk();
        i = i+1;
    end
    
    assert(~any(any(chunksO~=chunksO_)));
    assert(~any(any(chunks~=chunks_)));
    
    C.reset();
    
	chunksO_ = zeros(size(chunks));
    chunks_  = chunksO;
    i = 1;
    while C.hasNextChunk()
        [chunksO_(i,:) chunks_(i,:)] = C.getNextChunk();
        i = i+1;
    end
    
    assert(~any(any(chunksO~=chunksO_)));
    assert(~any(any(chunks~=chunks_)));        
    
    
    %% TEST 5
    chunksize = 1000;
    Len = 100;
    minchunksize = 10;
    overlap = 20;
    
    chunksO =[  1 100];
    chunks =[   1 100];
    
    C = mysort.util.Chunker(Len, 'chunkSize', chunksize, ...
                                 'chunkOverlap', overlap,...
                                 'minChunkSize', minchunksize);
	
	chunksO_ = zeros(size(chunks));
    chunks_  = chunksO;
    i = 1;
    while C.hasNextChunk()
        [chunksO_(i,:) chunks_(i,:)] = C.getNextChunk();
        i = i+1;
    end
    
    assert(~any(any(chunksO~=chunksO_)));
    assert(~any(any(chunks~=chunks_)));
    
    C.reset();
    
	chunksO_ = zeros(size(chunks));
    chunks_  = chunksO;
    i = 1;
    while C.hasNextChunk()
        [chunksO_(i,:) chunks_(i,:)] = C.getNextChunk();
        i = i+1;
    end
    
    assert(~any(any(chunksO~=chunksO_)));
    assert(~any(any(chunks~=chunks_)));            