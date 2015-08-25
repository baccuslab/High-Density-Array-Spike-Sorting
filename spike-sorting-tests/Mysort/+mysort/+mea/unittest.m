%% Definitions
% rootpath = '~/bel.svn/hima_internal/cmosmea_recordings/trunk';
% dpath = '/Roska/3Nov2011_DSGCs_rabbit';
% flist = {'Trace_id858_2011-11-03T17_16_27_11.stream.ntk'}; 
% fname_h5 = 'Trace_id858_2011-11-03T17_16_27_11.stream.h5';

rootpath = '/links/groups/hima/';
dpath = 'temporary/jamuelle/';
flist = {'Trace_id460_2010-04-15T10_21_21_34.stream.ntk'};
fname_h5 = 'Trace_id460_2010-04-15T10_21_21_34.stream.h5';

% image kan�le 130 (131?) zwei bits togglen neues frame
% matlab inkrementiert framenumber wenn bit togglet

load_chunk_size = 20000;
max_samples = 2000;
hpf = 500;
lpf = 3000;
%% Load NTK with NTK Loader
fullpath = fullfile(rootpath, dpath);

pathBuffer = pwd;
cd(fullpath)
ntk = initialize_ntkstruct(flist{1}, 'hpf', hpf, 'lpf', lpf);
[ntk2 ntk]=ntk_load(ntk, siz);
cd(pathBuffer);

% remove noisy channels and flat channels
ntk2=detect_valid_channels(ntk2,1);

% Plot some data
figure
subplot(1,4,1:3)
signalplotter(ntk2,'chidx',1:26, 'max_samples', max_samples);
subplot(1,4,4)
plot(ntk2.x, ntk2.y, '.')
axis equal
axis ij



%% Load H5 with CMOS MEA
ffileh5 = fullfile(rootpath, dpath, fname_h5);
MEA = mysort.mea.CMOSMEA(ffileh5);
ntk3 = MEA.initialize_ntkstruct('hpf', hpf, 'lpf', lpf);
[ntk4 ntk_in] = ntk_load(MEA, load_chunk_size);

%% COMPARISON between the two files

