function filtfilt_chunked(DS1, DS2, HD, filterChunkSize, prefilterScaleFactor)    
    T = size(DS1,1);
    assert(isa(DS1, 'mysort.ds.FilteredDataSourceInterface'), 'Must be a filterDS!');
%     channels = 1:size(DS1,2);

    chunker = mysort.util.Chunker(T, 'chunkSize', filterChunkSize,...
        'progressDisplay', 'console', 'chunkOverlap', 0, 'minChunkSize', 1000);
    
    tStart=tic;
    disp('Starting filtering...');    
    while chunker.hasNextChunk()
        chunk = chunker.getNextChunk();
        DS2(chunk(1):chunk(2),:) = int16(DS1(chunk(1):chunk(2),:)*prefilterScaleFactor);
    end
    tEnd = toc(tStart);
    fprintf('Total time for filtering: %.3f\n', tEnd);    
        
    
    
    %%%% OLD VERSION
%     fprintf('Starting forward prefiltering...\n');
%     
% %     if use_parfor
% %         disp('Initialising parallel pool...');
% %         matlabpool(4)
% %         disp('Starting filtering...');
% %         cstarts = chunker.chunks(:,1);
% %         cends = chunker.chunks(:,2);
% %         parfor i=1:size(chunker.chunks,1)
% %             chunk = [cstarts(i) cends(i)];
% %             cidx = cstarts(i):cends(i);
% %             R = func(chunk(1):chunk(2),channels);
% %             Y = filter(HD, R);
% %             Y = int16(Y*prefilterScaleFactor);
% %             DS2(cidx,:) = Y;
% %         end
% %     else
%         while chunker.hasNextChunk()
%             chunk = chunker.getNextChunk();
%             R = func(chunk(1):chunk(2),channels);
%             Y = filter(HD, R);
%             Y = int16(Y*prefilterScaleFactor);
%             DS2(chunk(1):chunk(2),:) = Y;
%         end
%         tEnd = toc(tStart);
%         fprintf('Total: %.3f\n', tEnd);
%         % init filter to reduce filter artifact at beginning
%         R = flipud(DS2(end-1000:end,:));   
%         Y = filter(HD, R);
% 
%         tStart=tic;
%         fprintf('Starting backward prefiltering...\n');
%         chunker = mysort.util.Chunker(T, 'chunkSize', filterChunkSize,...
%             'progressDisplay', 'console', 'chunkOverlap', 0, 'minChunkSize', 1000);
%         while chunker.hasNextChunk()
%             chunk = chunker.getNextChunk();
%             a = T-chunk(2)+1;
%             b = T-chunk(1)+1;
%             R = double(DS2(a:b,:))/prefilterScaleFactor;
%             R = flipud(R);
%             Y = filter(HD, R);
%             Y = flipud(Y);
%             Y = int16(Y*prefilterScaleFactor);
%             DS2(a:b,:) = Y;   % and overwrite forward filtered data
%         end          
% %     end
%     tEnd = toc(tStart);
%     fprintf('Total: %.3f\n', tEnd);            