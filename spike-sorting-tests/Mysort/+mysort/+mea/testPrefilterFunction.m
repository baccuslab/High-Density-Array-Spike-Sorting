% INIT
tmp_path = 'C:\LocalData\Michele\';
fname = 'Trace_id858_2011-11-03T17_16_27_11.stream.h5';
f_buffer_f = fullfile(tmp_path, [fname(1:end-2) '10_400_6000_1.h5']);

%%
P = struct();
P.hpf = 400;
P.lpf = 6000;
P.forder = 10;
P.chunkSize = 20000;
P.filtfilt = 1;
P.deflation = 1;
mysort.mea.prefilterH5Data(fullfile(tmp_path, fname), f_buffer_f, P);