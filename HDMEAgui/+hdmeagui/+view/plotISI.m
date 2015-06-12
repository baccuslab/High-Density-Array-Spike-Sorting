function plotISI(ax, GUI, t, C, GC)
    
    fprintf('  Plotting ISI...\n');
    temp = GUI.getCurrentTemplate();
    
    sr = GUI.WF.samplesPerSecond;
    tic
    if temp.getNExcludedWfs()>0
        tim = temp.getSourceSpikeTrain();
        mysort.plot.isi({tim}, 'figure', 0,...
            'axesHandles', ax, 'srate', sr,...
            'edges', 0:20:400, 'showIDs', 0, 'barParams', {'FaceColor', C.color_unselected_spikes});
    end
    tim = temp.getSpikeTrain();
    mysort.plot.isi({tim}, 'figure', 0,...
        'axesHandles', ax, 'srate', sr,...
        'edges', 0:20:400, 'showIDs', 0, 'barParams', {'FaceColor', C.color_selected_spikes});
    toc
    
    
    