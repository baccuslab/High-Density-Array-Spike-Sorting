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
loadings = [];
conversions = [];
for i=1:1 % length(cidx)
    c = cidx(i);
    comp = compressions(c);
    fname = sprintf(fstr, comp);

    ffile = fullfile(dpath, fname);
    fprintf('File: %s\n', fname);

    %%
    mea = mysort.mea.CMOSMEA(ffile);
    
    path = '/Sessions/Session0/sig';
    
    tic
    X = hdf5read(ffile, path);
    loadings(i) = toc;
    
%     tic
%     Xd = double(X);
%     conversions(i) = toc;
%     R = mea(1:2, 1:100);
%     disp(R);
    % mea(:)
end


%%
figure
subplot(2,1,1)
plot(compressions, loadings);
hold on

M=[
0               94.1347          25.021
1               26.4784          33.181
2               25.3418          36.062
3               24.4129          46.814
4               23.1222          42.338
5               22.9365          1*60+01.512
6               22.0552          1*60+44.759
7               21.7164          2*60+30.306
8               21.4802          8*60+18.336
9               21.3685         11*60+56.272];
plot(compressions, M(:,3), 'r')
legend('loading time', 'time to save')
ylabel('seconds')

subplot(2,1,2)
plot(compressions, M(:,2), 'g')
legend('compression ratio')
ylabel('ratio')
xlabel('compression number')