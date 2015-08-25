maxTau = 1;

nF = 4;
D = [ zeros(100, nF);
     [ (1:10)' (1:10)'*10 100*(1:10)' 1000*(1:10)'];
    zeros(100, nF); ];

nF = 2;
D = [ zeros(100, nF);
     [ (1:10)' (1:10)'*10];
    zeros(100, nF); ];

M = zeros(2*maxTau+1, nF, nF);
M(:,1,2) = maxTau:-1:-maxTau;
DC = mysort.sorters.DiscriminantFunctionContainer(D, M, maxTau);

figure;
plot(DC.DOvp)