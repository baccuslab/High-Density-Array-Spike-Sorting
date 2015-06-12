L = 1000000;

n = 1;
tstart = tic;
2+2;
t(n) = toc(tstart);
n = n+1; 

tstart = tic;
X = randn(L,100);
t(n) = toc(tstart);
n = n+1;

tstart = tic;
f = figure;
plot(X(:,1:2));
drawnow
close(f)
t(n) = toc(tstart);
n = n+1;

tstart = tic;
addpath('');
t(n) = toc(tstart);
n = n+1;

tstart = tic;
Y = filter(ones(100,1), 1, X);
t(n) = toc(tstart);
n = n+1;

tstart = tic;
Z = Y;
Z(end,end) = 2;
t(n) = toc(tstart);
n = n+1;

tstart = tic;
clear Y
t(n) = toc(tstart);
n = n+1;

tstart = tic;
Z(end+1,end) = 1;
t(n) = toc(tstart);
n = n+1;

tstart = tic;
clear Z
t(n) = toc(tstart);
n = n+1;

tstart = tic;
f = figure;
plot(X(:,1:2));
drawnow
close(f)
t(n) = toc(tstart);
n = n+1;


strs = {'Adding' 'Allocation' 'plotting 1', 'addpath' 'filter', 'copy', 'clear 1', 'elongation',  'clear 2', 'plotting 2'};
disp('###');
for i=1:length(strs)
    disp([strs{i} sprintf(': %.5fs', t(i))]);
end

dstr = datestr(now);
idx = strfind(dstr, ':');
dstr(idx) = '_';
save(dstr, 'strs', 't', 'dstr');
