function fh = sortAna(S, DS)
    nFh = 1; fh = [];
    notNoiseIdx = S.mergedMS.gdf(:,1)>1;
    fprintf('Spikes total: %d, NotNoiseSpikes: %d (%.4f%%)\n', size(S.mergedMS.gdf,1), sum(notNoiseIdx), 100*sum(notNoiseIdx)/size(S.mergedMS.gdf,1));
%     mysort.plot.clustering(S.fetXPCA(notNoiseIdx,:), S.mergedMS.gdf(notNoiseIdx,1));
    mysort.plot.clustering(S.fetXPCA, S.mergedMS.gdf);
    mysort.plot.figureName([ S.name ' MS Clustering']);
    fh(nFh) = gcf; nFh=nFh+1;
    
    mysort.plot.clusterPCA(S.spikesCut(notNoiseIdx,:), S.mergedMS.gdf(notNoiseIdx,1), S.C_time_cut);
    mysort.plot.figureName([ S.name ' MS Intra PCA']);
    fh(nFh) = gcf; nFh=nFh+1;

%     mysort.plot.waveforms(S.spikesCut(S.mergedMS.gdf(:,1)==1,:), 'nC', S.nC, 'linewidth', 3);
%     mysort.plot.figureName([ S.name ' MS Noise Spikes']);
%     fh(nFh) = gcf; nFh=nFh+1;
    
    mysort.plot.waveforms(S.mergedMS.templates(2:end,:), 'IDs', 2:size(S.mergedMS.templates,1), 'nC', S.nC, 'linewidth', 3);
    mysort.plot.figureName([ S.name ' MS Templates']);
    fh(nFh) = gcf; nFh=nFh+1;
    
%     mysort.plot.waveforms(S.mergedMS.templates(2:end,:), 'IDs', 2:size(S.mergedMS.templates,1), 'nC', S.nC, 'linewidth', 3, 'stacked', 0);
%     mysort.plot.figureName([ S.name ' MS Templates Single Plots']);
%     fh(nFh) = gcf; nFh=nFh+1;    
 
    mysort.plot.clusterProjection(S.spikesPW(notNoiseIdx,:), S.mergedMS.gdf(notNoiseIdx,1));
    mysort.plot.figureName([ S.name ' MS ClusterProj']);
    fh(nFh) = gcf; nFh=nFh+1;
     
    mysort.plot.isi(S.mergedMS.gdf, 'srate', S.srate, 'print2ms', 1, 'edges', (0:1:20)*S.srate/1000);
    mysort.plot.figureName([ S.name ' MS ISI']);
    fh(nFh) = gcf; nFh=nFh+1;
     
    mysort.plot.waveforms(S.spikesCut(notNoiseIdx,:), 'IDs', S.mergedMS.gdf(notNoiseIdx,1),...
        'nC', S.nC, 'stacked', 0);
    mysort.plot.figureName([ S.name ' MS Spike Wfs']);
    fh(nFh) = gcf; nFh=nFh+1;
     
    if nargin>1
        SpSoMS = mysort.spiketrain.SpikeSortingContainer('meanshift', S.mergedMS.gdf, ...
            'templateCutLeft', S.P.spike_cut_cutLeft, ...
            'templateCutLength', S.P.spike_cut_Tf, ...
            'templateWfs', mysort.wf.v2t(S.mergedMS.templates, S.nC),...
            'wfDataSource', DS, 'nMaxSpikesForTemplateCalc', 1000);
        DS.addSpikeSorting(SpSoMS); DS.setActiveSpikeSortingIdx('meanshift')
        DScell = {DS};
        chsp = 5000;
        if ~isempty(S.gdfbotm)
            DSbotm = mysort.ds.Matrix(DS, DS.samplesPerSecond, 'botm sorting');
            DScell = {DS, DSbotm};
            chsp = [5000 5000];
        end
        mysort.plot.SliderDataAxes(DScell, 'channelSpacers', chsp); 
        mysort.plot.figureName([ S.name ' MS Data w Templ']);
        fh(nFh) = gcf; nFh=nFh+1;
    end
end