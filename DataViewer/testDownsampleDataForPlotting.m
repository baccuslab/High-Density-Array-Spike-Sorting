%% ---- get Data 
DH = dataviewer.DataHandle('dbconfig', db.munk.DBconfig_Select);
X = DH.getData('trialIDs',1, 'tetrodeIDs', 3);

%% ---- calc
x = X(1,:);
x = x+50*sin(linspace(0,2*pi,length(x)));
[Y t] = downsample_for_plot(x, 4000);


%% ---- do plot
figure;

a(1) = subplot(2,1,1);
plot(t, Y, 'r');
hold on
plot(x, 'k');

a(2) = subplot(2,1,2);
plot(t, Y, 'k');

linkaxes(a, 'xy');