function handles = PlotTemplateButton(hObject, eventdata, handles)
    G = handles.GUI;
    T = G.T;
    GC = handles.GUI_CONFIG;
    C = handles.CONFIG;
    nC = T.wfs.nC; 
    Tf = GC.DisplaySpikesTf;
    
    fh = mysort.plot.figure('w', 600, 'h', 800');
    ax = axes();
    hold on
    
       %% TODO !!! MULTI ELECTRODE
    CL = T.wfs.dataSource.getChannelList();
    CL = CL(CL(:,2)==1,3:4)/1000;
    I = T.getIterator();
    while I.hasNext();
        t = I.next();
        tidx = I.idx;
        tid = t.id;
        [sTemp bTemp] = hdmeagui.view.templatesGetSpikesForPlotting(T, GC, tid);
        bTempM = mysort.wf.v2m(bTemp, nC);
        sTempM = mysort.wf.v2m(sTemp, nC);        
        bTempM = bTempM(:,1:end-1)/7;
        sTempM = sTempM(:,1:end-1)/7;

        if t.isAccepted()
            ls = '-';
            leg{tidx} = sprintf('T%d', tid);
        else
            ls = ':';
            leg{tidx} = sprintf('t%d', tid);
        end
        for i=1:nC
            idx = CL(i,1) + 15*(0:Tf-1)/Tf;
            y   = CL(i,2) + 1*sTempM(i,:);
            tempHandles(tidx) = plot(ax, idx, y, ls, 'color', mysort.plot.vectorColor(tid), 'linewidth', 2);
        end
    end
    legend(tempHandles, leg, 'location', 'NorthEastOutside');
    axis(ax, 'tight');
