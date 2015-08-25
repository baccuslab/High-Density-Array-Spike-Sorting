d = 2;
k = 5;
N = [150 117 227 137 134];
[X P T] = mysort.clustering.generateTestData('dim', nDims, 'k', k, 'mu_spread', 50, 'N', N);


% n = 400; k = 7; d = 2; c = 3; e = 10;
% fprintf('Sampling a %d-component %d-dimensional %d-separated mixture....\n',k,d,c);
% [X,T] = greedyEM.mixgen(n,n,k,d,c,e);
fprintf('------------------\n');

fprintf('Running 10 times normal EM with k-means initialization\n');
normal_em=[];
for i=1:10
  [W,M,R,Tlogl] = greedyEM.em(X,T,k,0,0,0);
  normal_em = [normal_em Tlogl];
end
MEM = M;
fprintf('Average log-likelihood %f with std. dev. %f best run: %f \n',mean(normal_em),std(normal_em), max(normal_em));

     max_k = 9;
fprintf('Running greedy EM\n');
[W,M,R,Tlogl] = greedyEM.em(X,T,max_k,10,1,0);
title('Mixture model');


figure(2); clf;plot(Tlogl,'-o');hold on;
xlabel 'number of components'
ylabel 'log-likelihood of test set'
plot(repmat(mean(normal_em),max_k,1),'r')
plot(repmat(mean(normal_em)+std(normal_em),max_k,1),'r--')
plot(repmat(max(normal_em),max_k,1),'g')
plot(repmat(mean(normal_em)-std(normal_em),max_k,1),'r--')

[Tlogl,best_k] = max(Tlogl);
fprintf('Best number of components according to cross-validation: %d (yielding log-likelihood %f ) \n',best_k,Tlogl);


legend('Greedy EM','mean of normal EM', 'one standard deviation margins of normal EM','best result of normal EM',4);
title('Log-likelihood plots');
figure(1)
hold on
plot(P.mu(:,1), P.mu(:,2), 'rx', 'markersize', 20, 'linewidth', 4);
plot(M(1:best_k,1), M(1:best_k,2), 'bo', 'markersize', 20, 'linewidth', 4);
plot(MEM(:,1), MEM(:,2), 'kd', 'markersize', 20, 'linewidth', 4);


