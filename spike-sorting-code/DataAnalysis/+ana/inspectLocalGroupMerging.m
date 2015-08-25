function inspectLocalGroupMerging(S, local_nC)
    mysort.plot.waveformsVertical(S.ALI.spikeAligned.wfs, 'stacked', 0, 'IDs', S.CLU.clustering.ids, 'nC', local_nC)
    mysort.plot.figureTitle('Aligned spikes, cluster IDs');
    
    idx = S.ALI.spikeAligned.alignIdx;
    mysort.plot.waveformsVertical(S.ALI.spikeAligned.wfs, 'stacked', 0, 'IDs', S.MER.clusteringMerged.ids(idx), 'nC', local_nC)
    mysort.plot.figureTitle('Aligned spikes, merged IDs');
    
