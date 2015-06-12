%% Init
defs = mysort.mea.test.definitions();
h5f  = [defs.testfile_location defs.testfile_h5];
ntkf =  [defs.testfile_location defs.testfile_ntk];

TT = [10^5 10^6 2*10^6];
reps = 5;
TIMES = zeros(2, length(TT), reps);
for t = 1:length(TT)
    T = TT(t);
    for rep = 1:reps
        %% INIT NTk
        pathBuffer = pwd;
        cd(defs.testfile_location);
        ntk = initialize_ntkstruct(defs.testfile_ntk, 'hpf', 500, 'lpf', 3000);


        %% INIT Mea
        mea = mysort.mea.CMOSMEA(h5f, 'hpf', 500, 'lpf', 3000);


        %% Load NTk
        fprintf('Loading NTK...')
        tic
        [ntk2 ntk] = ntk_load(ntk, T);
        cd(pathBuffer);
        TIMES(1, t, rep) = toc;

        %% Load H5
        fprintf('Loading H5...')
        tic
        Xmea = mea(1:T, 1:92);
        TIMES(2, t, rep) = toc;
    end
end

%save('time_run2', 'TIMES', 'rep', 'TT', 'defs')

%%
figure;
ah(1) = subplot(1,2,1);
boxplot(squeeze(TIMES(1,:,:))')
title(sprintf('Loading times NTK, N=%d', rep));
set(gca, 'xtick', 1:length(TT), 'xticklabel', TT)
xlabel('Samples to load')
ylabel('Time to load (sec)');

ah(2) = subplot(1,2,2);
boxplot(squeeze(TIMES(2,:,:))', 'colors', 'g')
title(sprintf('Loading times CMOSMEA, N=%d', rep))
set(gca, 'xtick', 1:length(TT), 'xticklabel', TT)
xlabel('Samples to load')
linkaxes(ah, 'xy');
% mysort.plot.savefig(gcf, 'time_run1')
