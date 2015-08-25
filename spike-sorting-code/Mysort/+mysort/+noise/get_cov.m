

%% Compute noise epochs common with n nearest electrodes
NEN = NE;
for ch = 1:nC
    % find k nearest electrodes
    neighbors = mysort.mea.nearestElectrodes(channel_x, channel_y, ch, 6);
    for n = 1:length(neighbors)
        NEN{ch} = mysort.epoch.intersect(NEN{ch}, NE{neighbors(n)});
    end
end


%%
nXc = (nC*(nC-1))/2;
D = zeros(1,nXc);count = 1;
Tf = 30; maxLag = Tf-1;
autocovs = zeros(nC, 2*Tf-1);
autocov_norms = zeros(nC,1);
xcovs = zeros(nXc, 2*Tf-1);
IDX = zeros(nXc,2);


for ch1 = 1:nC
    p1 = [channel_x(ch1) channel_y(ch1)];
    [ac lc] = mysort.util.xcorr_in_epochs(Y(ch1,:), NE{ch1}, maxLag, maxLag);
    autocovs(ch1,:) = ac;
    autocov_norms(ch1) = lc;
    
    for ch2 = ch1+1:nC 
        p2 = [channel_x(ch2) channel_y(ch2)];
        d = norm(p1-p2);
        IDX(count,:) = [ch1 ch2];
        D(count) = d; 
        if d < 40
            cNE = mysort.epoch.intersect(NE{ch1},NE{ch2});
            cNE = mysort.epoch.removeShort(cNE, 10*Tf);
            cNE = cNE(1:3:end,:);
            [xc lc] = mysort.util.xcorr_in_epochs(Y([ch1 ch2],:), cNE, maxLag, maxLag);
        else
            xc = zeros(1, 2*Tf-1);
        end
        xcovs(count,:) = xc;   
        count=count+1;
        disp(count)
    end
end
% %%
% save('C:\LocalData\Michele\marching_square_buffer_covs', 'D', 'Tf', 'autocovs', ...
%      'autocov_norms', 'D', 'xcovs', 'IDX', 'smad', 'NE');
% %%
% noiseepochs1 = mysort.epoch.flip(spikeepochs1, size(X,2));
% [smad2, spikeepochs2] = ana.douglas.estimateSigma(X(k2,:), Tf, thr1);
% noiseepochs2 = mysort.epoch.flip(spikeepochs2, size(X,2));
% 
% spikeepochs1_idx = mysort.epoch.toIdx(spikeepochs1);
% spikeepochs2_idx = mysort.epoch.toIdx(spikeepochs2);
% noiseepochs1_idx = mysort.epoch.toIdx(noiseepochs1);
% noiseepochs2_idx = mysort.epoch.toIdx(noiseepochs2);
% commonnoiseepochs = mysort.epoch.intersect(noiseepochs1,noiseepochs2);
% commonnoiseepochs_idx = mysort.epoch.toIdx(commonnoiseepochs);
% unionspikeepochs  = mysort.epoch.merge( [spikeepochs1; spikeepochs2]);
% unionspikeepochs_idx = mysort.epoch.toIdx(unionspikeepochs);