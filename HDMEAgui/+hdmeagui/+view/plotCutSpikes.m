function plotCutSpikes(ax, T, t, C, GC)
   
    % prepare spikes for plotting
    tic
    fprintf('    Preparing Spikes for plotting...\n');
    [sTemp bTemp selSpikes unselSpikes uTemp] = hdmeagui.view.templatesGetSpikesForPlotting(T, GC, t, GC.NPlotChans);
    toc 
    
    if isempty(unselSpikes)
        fprintf('  Plotting Spikes...\n');
        tic;
        set(ax, 'NextPlot', 'replacechildren');
        plot(ax, selSpikes', 'color', C.color_selected_spikes, 'linewidth',1);
        set(ax, 'NextPlot', 'add');
%         plot(ax, bTemp, 'color', 'k', 'linewidth', 2);
        plot(ax, sTemp, 'color', C.color_selected_template, 'linewidth', 2);
        axis(ax, 'tight')    
        toc 
    else
        fprintf('Plotting...\n');
        % plot spikes
        tic
        set(ax, 'NextPlot', 'replacechildren');
        plot(ax, unselSpikes', 'color', C.color_unselected_spikes, 'linewidth',1);
        set(ax, 'NextPlot', 'add');
        plot(ax, selSpikes', 'color', C.color_selected_spikes, 'linewidth',1);        
        plot(ax, bTemp, C.color_unselected_template, 'linewidth', 3);
        plot(ax, sTemp, 'color', C.color_selected_template, 'linewidth', 2);
        plot(ax, uTemp, 'color', 'b', 'linewidth', 2);
        axis(ax, 'tight')            
        toc
    end
    xlabel(ax, 'time in every channel')
    ylabel(ax, 'amplitude [std noise]');
    

end