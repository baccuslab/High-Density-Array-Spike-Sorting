function handles = PlotLoadRewardReaction(hObject, eventdata, handles)
    P = dataviewer.util.getHandleValues(handles);
    srate = handles.DH.getSamplesPerSecond();
    Q = handles.DH.query(['SELECT t.trialidx, t.ok, t.reward_time, t.response_time, t.load, t.length '...
                          'FROM trial as t WHERE t.id IN (' util.dlmstring(P.trialIDs) ') '...
                          'ORDER BY t.trialidx ASC']);
	% Get every column individually since only null columns will be omitted
	% by cell2mat
    T = -1*ones(size(Q));
	T(:,1:2) = cell2mat(Q(:,1:2));
    for i=1:size(T,1)
        for j=3:6
            if ~isempty(Q{i,j})
                T(i,j) = Q{i,j};
            end
        end
    end
    
    ok  = T(:,2) == 1;
    rew = T(:,3) > -1;
    
	fig = mysort.plot.figure('name', 'TrialDetailPlot');
    mysort.plot.figureTitle(P.figureTitle);
    ah = mysort.plot.subplot([2,1]);
    
    l0 = T(:,5) == -1;
    l1 = T(:,5) == 1;
    l2 = T(:,5) == 2;
    l3 = T(:,5) == 3;
    l4 = T(:,5) == 4;
    
    axes(ah(1)); hold on;
    legendstr = {};
    if any(l0)
        plot(T(l0,1), 0*ones(sum(l0), 1), 'kx');
        legendstr = {'Defect'};
    end
    
    colors = {'b.', 'bx', 'r.', 'rx', 'g.', 'gx', 'c.', 'cx'};
    for L = 1:4
        l = T(:,5) == L;
        if any(l &  rew)
            plot(T(l &  rew,1), L*ones(sum(l &  rew), 1), colors{2*L-1});
            legendstr = [legendstr ['Load ' num2str(L) ' (rewarded)']];
        end
        if any(l & ~rew)
            plot(T(l & ~rew,1), L*ones(sum(l & ~rew), 1), colors{2*L  });
            legendstr = [legendstr ['Load ' num2str(L) ' (not rew.)']];
        end    
    end
    set(ah(1), 'xticklabel', []);
    axis tight
    set(ah(1), 'ylim', [-.5 4.5]);

    ylabel('Load');
    legend(legendstr);
    
    axes(ah(2)); hold on
    plot(T( ok,1),    T(ok,6)/srate, 'b.');
    plot(T(~ok,1),    zeros(sum(~ok),1), 'bx');    
    plot(T( rew,1),   T( rew,4)/srate, 'r.');
    plot(T(~rew,1),   T(~rew,4)/srate, 'rx');
    legend('length (trial ok)', 'length (trial not ok)', 'reaction time (rewarded)', 'reaction time (not rew.)', 'location', 'East');
    ylabel('time [s]');
    xlabel('trial idx');
    axis tight
    linkaxes(ah, 'x');
    
    
    handles.new_figure_handles = fig;