function compareSortingSlider(DS1, DS2, DS3, DS4, S)
    SpSoMS = mysort.spiketrain.SpikeSortingContainer('meanshift', S.clusteringMatched.gdf, ...
            'templateCutLeft', S.P.spikeCutting.cutLeft, ...
            'templateCutLength', S.P.spikeCutting.Tf, ...
            'templateWfs', mysort.wf.v2t(S.clusteringMatched.templates, S.nC),...
            'wfDataSource', DS1, 'nMaxSpikesForTemplateCalc', 1000);
    DS1.addSpikeSorting(SpSoMS); DS1.setActiveSpikeSortingIdx('meanshift')
    A1 = ana.moritzheimdahl.ArtefactDetector(DS1, S.P.artefactDetection.width, S.P.artefactDetection.threshold);
    pSort  = [1 0];
    DScell = {DS1 A1};
    noisestd = sqrt(mean(diag(S.noise.C_time_cut)));
    chsp = 3*noisestd*[1 0];
    if ~isempty(S.botm.gdf)
        % Botm
        SpSoBOTM = mysort.spiketrain.SpikeSortingContainer('botm', S.botm.gdf, ...
            'templateCutLeft', S.P.spike_cut_cutLeft, ...
            'templateCutLength', S.P.spike_cut_Tf, ...
            'templateWfs', mysort.wf.v2t(S.templatesAfterBOTM, S.nC),...
            'wfDataSource', DS2, 'nMaxSpikesForTemplateCalc', 1000);
        DS2.addSpikeSorting(SpSoBOTM); DS2.setActiveSpikeSortingIdx('botm')            
        % Residual
        DS3.addSpikeSorting(SpSoBOTM); DS3.setActiveSpikeSortingIdx('botm') 
        DS3.bReturnSortingResiduals = 1;
        A2 = ana.moritzheimdahl.ResidualArtefactDetector(DS3);
        
        % botm cleaned
        SpSoBOTMCl = mysort.spiketrain.SpikeSortingContainer('botmcl', S.gdfbotmcleaned, ...
            'templateCutLeft', S.P.spike_cut_cutLeft, ...
            'templateCutLength', S.P.spike_cut_Tf, ...
            'templateWfs', mysort.wf.v2t(S.templatesAfterBOTMcleaned, S.nC),...
            'wfDataSource', DS4, 'nMaxSpikesForTemplateCalc', 1000);
        DS4.addSpikeSorting(SpSoBOTMCl); DS4.setActiveSpikeSortingIdx('botmcl')            
        
        pSort  = [1    0   1     0   0   1];
        DScell = {DS1, A1, DS2, DS3, A2, DS4};
        chsp = 6*noisestd*[1 0 1 1 0 1];
    end
    
    mysort.plot.SliderDataAxes(DScell, 'channelSpacers', chsp, 'plotSortings', pSort, 'timeIn', 'samples'); 
    mysort.plot.figureName([ S.name ' MS Data w Templ']);