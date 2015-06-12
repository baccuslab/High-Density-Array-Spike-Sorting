%%
x = [0:.1:8*pi];
y = [2*sin(x);
     2*cos(x)];
sf = mysort.wf.mSincfun(x, y);
xi1 = 0:.05:2*pi;
yi1 = sf(xi1);

xi2 = 0:.2:2*pi;
yi2 = sf(xi2);

xi3 = 8+(0:.2:4*pi);
yi3 = sf(xi3);

figure;
subplot(2,1,1)
plot(x, y(1,:), '.-b');
hold on
plot(xi1, yi1(1,:), '.-r');
plot(xi2, yi2(1,:), '.-c');
plot(xi3, yi3(1,:), '.-g');
subplot(2,1,2)
plot(x, y(2,:), '.-b');
hold on
plot(xi1, yi1(2,:), '.-r');
plot(xi2, yi2(2,:), '.-c');
plot(xi3, yi3(2,:), '.-g');


%%
x = [0:.1:8*pi];
y = [x;
     -x];
sf = mysort.wf.mSincfun(x, y);
xi1 = 0:.05:2*pi;
yi1 = sf(xi1);

xi2 = 0:.2:2*pi;
yi2 = sf(xi2);

xi3 = 8+(0:.2:4*pi);
yi3 = sf(xi3);

figure;
subplot(2,1,1)
plot(x, y(1,:), '.-b');
hold on
plot(xi1, yi1(1,:), '.-r');
plot(xi2, yi2(1,:), '.-c');
plot(xi3, yi3(1,:), '.-g');
subplot(2,1,2)
plot(x, y(2,:), '.-b');
hold on
plot(xi1, yi1(2,:), '.-r');
plot(xi2, yi2(2,:), '.-c');
plot(xi3, yi3(2,:), '.-g');

%%
figure
x = -10:.05:10;
y = mysort.util.sinc0(x);
plot(x,y)
