function S = loadLocalGroupSorting(spath, rname, groupidx)
    S.rname = rname;
    S.groupidx = groupidx;
    S.pref = [spath '\group' sprintf('%03d', groupidx) '\' rname];
    S.CUT = load([S.pref '.040spikes_cut.mat']);
    S.ALI = load([S.pref '.050spikes_aligned.mat']);
    S.PW = load([S.pref '.070spikes_prewhitened.mat']);
    S.FET = load([S.pref '.080spikes_features.mat']);
    S.CLU = load([S.pref '.090clusters_meanshift.mat']);
    S.TM = load([S.pref '.100botm_matching.mat']);
    S.MER = load([S.pref '.110clusters_meanshift_merged.mat']);
   
    
    
