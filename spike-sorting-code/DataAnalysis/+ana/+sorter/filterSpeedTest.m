%         if 0
%             tic
%             Xoldhd = filtfilthd(hd, double(X));
%             toc
%         
%             tic
%             Xmat1 = zeros(size(X));
%             for c=1:size(X,2)
%                 Xmat1(:,c) = conv2(double(X(:,c)), mysort.mea.filter_gpu()','same');
%             end
%             toc
%             
%             tic
%             gpuF = gpuArray(mysort.mea.filter_gpu());
%             X2 = zeros(size(X));
%             for c=1:size(X,2)
%                 X2(:,c) = gather(conv2(gpuArray(double(X(:,c))), gpuF','same'));
%             end
%             toc
%             
%             
%             tic
%             X3 = conv2(double(X), mysort.mea.filter_gpu()','same');
%             toc
%             
%             tic
%             gpuF = gpuArray(mysort.mea.filter_gpu());
%             X4 = gather(conv2(gpuArray(double(X)), gpuF','same'));
%             toc    
%             
%             figure; plot(X(1:100,:),'k');
%             hold on; plot(Xoldhd(1:100,:), 'b'); 
%             hold on; plot(Xmat1(1:100,:), 'm');  
%             hold on; plot(X2(1:100,:), 'g');  
%             hold on; plot(X3(1:100,:), 'r');    
%         end