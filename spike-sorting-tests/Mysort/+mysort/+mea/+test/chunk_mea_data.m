if 1
    clear mea
    clear mysort.mea.CMOSMEA
end

fstr = 'Trace_id843_2011-12-06T11_27_53_5.stream.h5';
mea = mysort.mea.CMOSMEA(fstr);



[nT nC] = size(mea);
chunkSize = 500000;

% For test purposes
nT = 1000000;
chunkSize = 200000;

chunker = mysort.util.Chunker(nT, 'overlap', 0, ...
                    'chunkSize', chunkSize,...
                    'progressDisplay', 'console');
tic                
while chunker.hasNextChunk()
    [chunk_overlap chunk chunkLen] = chunker.getNextChunk();
    X = mea(chunk, 1:nC);
end
toc

