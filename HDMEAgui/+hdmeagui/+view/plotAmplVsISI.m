function plotAmplVsISI(ax, GUI, t, C, GC)    
    maxT = 150*GUI.WF.samplesPerSecond/1000;
    set(ax, 'Nextplot', 'replace');
        
    a = GUI.getCurrentTemplateSourceProjections();
    temp = GUI.getCurrentTemplate();
    tim = temp.getSourceSpikeTrain();
    ppar = {'line', 'none', 'marker', 'x', 'color', C.color_unselected_spikes};
    mysort.plot.amplitudeVsIsi(a,tim, 'axesHandles', ax, 'maxIsi', maxT, ...
        'plotParas', ppar, 'srate', GUI.WF.samplesPerSecond, 'd3', 0)
    
    if temp.getNExcludedWfs()>0
        set(ax, 'Nextplot', 'add');
        a = GUI.getCurrentTemplateProjections();
        tim = temp.getSpikeTrain();
        ppar = {'line', 'none', 'marker', 'x', 'color', C.color_selected_spikes};
        mysort.plot.amplitudeVsIsi(a,tim, 'axesHandles', ax, ...
            'plotParas', ppar, 'maxIsi', maxT, 'srate', GUI.WF.samplesPerSecond, 'd3', 0)
    end
    
    xlabel(ax, 'isi [ms]');
    ylabel(ax, 'projection');
    set(ax, 'box', 'off');
    
    