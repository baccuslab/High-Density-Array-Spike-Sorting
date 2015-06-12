function handles = PlotRates(hObject, eventdata, handles)
    P = dataviewer.util.getHandleValues(handles);
    srate = handles.DH.getSamplesPerSecond();
    Q = handles.DH.query([...
'SELECT t.trialidx, u.id, u.algoid, (cast(count(*) AS float))/(cast(t.length AS float)) FROM trial as t '...
'JOIN spike as s ON (s.trial = t.id) '...
'JOIN unit as u ON (u.id = s.unit) '...
'WHERE t.id IN (' util.dlmstring(P.trialIDs) ') '...
'AND u.id IN (' util.dlmstring(P.unitIDs) ') '...
'GROUP BY t.trialidx, u.id, u.algoid, t.length '...
'ORDER BY t.trialidx, u.algoid']);
	Q = cell2mat(Q);
    
    T = unique(Q(:,1));
    U = unique(Q(:,2));
    A = unique(Q(:,3));
    R = zeros(length(T), length(U));    
    for i=1:size(Q,1)
        u = find(Q(i,2)==U,1);
        t = find(Q(i,1)==T,1);
        R(t, u) = Q(i,4);
    end
    
    R = R*srate;
    
	fig = mysort.plot.figure('name', 'RatePlot');
    mysort.plot.figureTitle(P.figureTitle);
    ah = mysort.plot.subplot([1,1]);
   
    
    axes(ah(1)); hold on;
    
    legendstr = {};
    for u=1:length(U)
        plot(T, R(:,u), 'color', mysort.plot.vectorColor(A(u)));
        legendstr = [legendstr {num2str(U(u))}];
    end
    legend(legendstr);

    ylabel('rate [Hz]');
    xlabel('trial idx');
    axis tight
    
    
    handles.new_figure_handles = fig;