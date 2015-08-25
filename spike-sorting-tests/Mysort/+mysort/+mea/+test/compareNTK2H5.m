%% Init
defs = mysort.mea.test.definitions();
h5f  = [defs.testfile_location defs.testfile_h5];
ntkf =  [defs.testfile_location defs.testfile_ntk];


%% Load NTk
siz=40000;
pathBuffer = pwd;
cd(defs.testfile_location);
ntk = initialize_ntkstruct(defs.testfile_ntk, 'hpf', 500, 'lpf', 3000);
[ntk2 ntk]               = ntk_load(ntk, siz);
cd(pathBuffer);

%% Load H5
mea = mysort.mea.CMOSMEA(h5f, 'hpf', 500, 'lpf', 3000);

%% Compare Raw
Xmea_raw = mea.getRawData(1:100, 1:10, 1, 'channels', 'all')';
Xntk_raw = ntk.data(1:10, 1:100);
spacer = mysort.plot.mc(double(Xmea_raw), 'color', {'k'}, 'linewidth', 2);
hold on
mysort.plot.mc(double(Xntk_raw), 'figure', 0, 'color', {'r'}, 'spacer', spacer);

%% Compare Filtered
size(ntk2.sig)
size(mea)
Xmea = mea(1:10000, 1:10);
Xntk = ntk2.sig(1:10000,1:10);
% sum(sum(abs(Xmea-Xntk)))

spacer = mysort.plot.mc(Xmea', 'color', {'k'});
hold on
mysort.plot.mc(Xntk', 'figure', 0, 'color', {'r'}, 'spacer', spacer);

figure; hist(Xmea(:)./Xntk(:),100)