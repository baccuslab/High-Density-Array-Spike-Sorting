 function fh = sortAnaBotm(S, DS)
    nFh = 1; fh = [];
%     notNoiseIdx = S.gdfms(:,1)>1;
%     mysort.plot.clustering(S.fetXPCA(notNoiseIdx,:), S.gdfms(notNoiseIdx,1));
%     mysort.plot.figureName([ S.name ' MS Clustering']);
%     fh(nFh) = gcf; nFh=nFh+1;
    
%     mysort.plot.clusterPCA(S.spikesCut(notNoiseIdx,:), S.gdfms(notNoiseIdx,1), S.C_time_cut);
%     mysort.plot.figureName([ S.name ' MS Intra PCA']);
%     fh(nFh) = gcf; nFh=nFh+1;
%     
    mysort.plot.waveforms(S.templates4BOTM, 'IDs', 1:size(S.templates4BOTM,1), 'nC', S.nC, 'linewidth', 3);
    mysort.plot.figureName([ S.name ' BOTM used Templates']);
    fh(nFh) = gcf; nFh=nFh+1;
    mysort.plot.waveforms(S.templates4BOTM, 'IDs', 1:size(S.templates4BOTM,1), 'nC', S.nC, 'linewidth', 3, 'stacked', 0);
    mysort.plot.figureName([ S.name ' BOTM used Templates']);
    fh(nFh) = gcf; nFh=nFh+1;    
%  
%     mysort.plot.clusterProjection(S.spikesPW(notNoiseIdx,:), S.gdfms(notNoiseIdx,1));
%     mysort.plot.figureName([ S.name ' MS ClusterProj']);
%     fh(nFh) = gcf; nFh=nFh+1;
     
    mysort.plot.isi(S.gdfbotmcleaned, 'srate', S.srate, 'print2ms', 1);
    mysort.plot.figureName([ S.name ' BOTM ISI']);
    fh(nFh) = gcf; nFh=nFh+1;
     
%     mysort.plot.waveforms(S.spikesCut(notNoiseIdx,:), 'IDs', S.gdfms(notNoiseIdx,1),...
%         'nC', S.nC, 'stacked', 0);
%     mysort.plot.figureName([ S.name ' MS Spike Wfs']);
%     fh(nFh) = gcf; nFh=nFh+1;