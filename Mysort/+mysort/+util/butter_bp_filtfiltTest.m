srate = 20000;
order = 10;
L = 1000;
x  = mysort.util.generateSinus(6000, 90, srate, L);
X = [x' -2*x'];

Y = mysort.util.butter_bp_filtfilt(X, order, 300, 8000, srate);

%%
figure;

ah = subplot(2,1,1);
plot(X(:,1), '.-b');
hold on
plot(Y(:,1), '--r');

ah(2) = subplot(2,1,2);
plot(X(:,2), '.-b');
hold on
plot(Y(:,2), '--r');
linkaxes(ah, 'xy')
%%
srate = 20000;
order = 10;
L = 1000;
x  = 1 + .7*mysort.util.generateSinus( 50,  0, srate, L) ...
       +    mysort.util.generateSinus(120, 90, srate, L) ...
       +    mysort.util.generateSinus(505, 90, srate, L) ...
       +    mysort.util.generateSinus(512, 90, srate, L) ...
       +    mysort.util.generateSinus(6000, 90, srate, L) ...
       +    mysort.util.generateSinus(8000, 90, srate, L) ...
       +    mysort.util.generateSinus(9000, 90, srate, L);
X = [x' -2*x'];

burnin = 200;
Y = mysort.util.butter_bp_filtfilt(X, order, 300, 8000, srate);

%%
figure;

ah = subplot(2,1,1);
plot(X(:,1), '.-b');
hold on
plot(Y(:,1), '--r');

ah(2) = subplot(2,1,2);
plot(X(:,2), '.-b');
hold on
plot(Y(:,2), '--r');
linkaxes(ah, 'xy')