%%
fname = 'Trace_id843_2011-12-06T11_27_53_5.stream.h5';
plist = 'H5P_DEFAULT';
rmode = 'H5F_ACC_RDONLY';
channelSet = [1:92];
nC = length(channelSet);
T = 100000;

%% Read from mea object
tic
mea = mysort.mea.CMOSMEA(fstr);
X0 = mea(1:T, channelSet);
toc