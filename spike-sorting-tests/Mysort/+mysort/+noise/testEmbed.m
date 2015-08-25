X = randn(3,100);

XX = [X; [X(:,2:end) zeros(3,1)]; [X(:,3:end) zeros(3,2)]]; 

C = XX*XX';

figure;
imagesc(C)