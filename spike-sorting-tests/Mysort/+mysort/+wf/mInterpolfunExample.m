% % %%
x = 0:.1:8*pi;
y = 2*sin(x);
sf = mysort.wf.mInterpolfun(x, y);
xi1 = 0:.05:2*pi;
yi1 = sf(xi1);

xi2 = 0:.2:2*pi;
yi2 = sf(xi2);

xi3 = 8+(0:.2:4*pi);
yi3 = sf(xi3);

figure;
plot(x, y, '.-b');
hold on
plot(xi1, yi1, '.-r');
plot(xi2, yi2, '.-c');
plot(xi3, yi3, '.-g');

%%
nTest = 4000;
tic
for i=1:nTest
    y = sf(pi);
end
t = toc;
fprintf('Time per interpolation: %f\n', t/nTest);


%%
% y = sin(0:.01:2*pi);
% figure; plot(y)
% y(1,1000) = 0;
% y(313:930) = .3*sin(linspace(pi, 2*pi, 930-313+1));
% figure; plot(y)
% y(100:220) = y(100:220) + sin(linspace(0, pi, 220-100+1));
% figure; plot(y)
% y(159-3:159+3) = y(159-3:159+3).*[1.05 1.1 1.2 1.35 1.2 1.1 1.05]; 
% figure; plot(y)