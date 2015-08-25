%%
if 0
    clear mea
    clear mysort.mea.CMOSMEA
end

%% Init
dpath = 'C:\LocalData\H5files';
fstr = 'Trace_id460_2010-04-15T10_21_21_34_compression_%d.stream.h5';

compressions = 0:9;

cidx = 1:length(compressions);

c = cidx(1);
comp = compressions(c);
fname = sprintf(fstr, comp);

ffile = fullfile(dpath, fname);
fprintf('File: %s\n', fname);

%%
mea = mysort.mea.CMOSMEA(ffile);
t = 100;
l = 30;
mea(5, t-30:t+30)
CL = mea.getChannelList();
[nC b] = size(mea)
mea.getChannelNr()'

%%
% nC = 2;
figure;
subplot(3,1,1)
plot(mea.getRawData('channels', 1:nC, 't1', 1, 't2', 100))
subplot(3,1,2)
plot(mea(1:nC, 1:100))
subplot(3,1,3)
% hpf: 500 lpf: 3000
% in hidense_filter.m line 52 (that is line 75 after davids last commit) the first bbp
