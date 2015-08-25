
function ah = alignment(R, varargin)
    defs = mysort.util.defs();
    P.figureHandle = [];
    P.doNotLabel = [defs.tp defs.tpo];
    P.mark = [defs.fp defs.fn defs.fno defs.cl defs.clo];
    P.connectors = true;
    P.axesHandle = [];
    P.stackGT = true;
    P.restrict2Time = []; % two element array, [startTime endTime]
    P = mysort.util.parseInputs(P, varargin, 'error');
    
    if isempty(P.figureHandle) && isempty(P.axesHandle)
        P.figureHandle = mysort.plot.figure('width', 1200, 'height', 600);
    end
    if isempty(P.axesHandle)
        P.axesHandle = axes(); 
    else
        axes(P.axesHandle)
    end    
    hold on
    
    if P.stackGT
        plotStackedGT();
    else
        plotStackedAssociations();
    end
    ylabel('Units in Sorted Spikes (St2)    vs.     Units in Ground Truth (St1)');
    xlabel('Time');
    ah = P.axesHandle;
    %----------------------------------------------------------------------
    function plotStackedAssociations()
        error('not implemented yet');
        if P.connectors
            for i = 1:size(R.ALI,1)
                for j = 1:size(R.ALI,2)
                    for k = 1:size(R.ALI{i,j}, 1)
                        ab = R.ALI{i,j}(k,:);
                        a = R.St1{i}(ab(1));
                        b = R.St2{j}(ab(2));
                        plot([a b], [i-.5 -i+.5], 'k:')
                    end
                end
            end
        end

%         for st1=1:length(R.St1)
%             plot(R.St1{st1}', st1*ones(1, length(R.St1{st1}))-.5, '.', 'color', mysort.plot.vectorColor(st1), 'markersize', 14);
%             for i=1:length(R.St1{st1})
%                 if any(R.spikeLabel1{st1}(i) == P.mark)
%                     plot(R.St1{st1}(i), st1-.5, 'rd', 'markersize', 10, 'linewidth', 2);
%                 end
%                 if ~any(R.spikeLabel1{st1}(i) == P.doNotLabel)
%                     text(R.St1{st1}(i), st1-.2, defs.labelID2String{R.spikeLabel1{st1}(i)});
%                 end
%             end
%         end
%         for st2=1:length(R.St2)
%             plot(R.St2{st2}', -st2*ones(1, length(R.St2{st2}))+.5, '.', 'color', mysort.plot.vectorColor(st2+length(R.St1)), 'markersize', 14);
%             for i=1:length(R.St2{st2})
%                 if any(R.spikeLabel2{st2}(i) == P.mark)
%                     plot(R.St2{st2}(i),-st2+.5, 'rd', 'markersize', 10, 'linewidth', 2);
%                 end                        
%                 if ~any(R.spikeLabel2{st2}(i) == P.doNotLabel)
%                     text(R.St2{st2}(i),-st2+.2, defs.labelID2String{R.spikeLabel2{st2}(i)});
%                 end
%             end        
%         end    

        axis tight
        xlim = get(gca, 'xlim');
        plot(xlim +[-10 10], [0 0], 'k-', 'linewidth', 2);
        axis tight
        set(gca, 'ylim', [-length(R.St2) max(1, length(R.St1))]); % make sure this is not [0 0] !
    end
    
    %----------------------------------------------------------------------
    function plotStackedGT()
        if P.connectors
            for i = 1:size(R.ALI,1)
                for j = 1:size(R.ALI,2)
                    for k = 1:size(R.ALI{i,j}, 1)
                        ab = R.ALI{i,j}(k,:);
                        a = R.St1{i}(ab(1));
                        b = R.St2{j}(ab(2));
                        if ~isempty(P.restrict2Time)
                            idx = a>P.restrict2Time(1) & a<P.restrict2Time(2) & b>P.restrict2Time(1) & b<P.restrict2Time(2);
                            a = a(idx);
                            b = b(idx);
                        end
                        if ~isempty(a)
                            plot([a b], [i-.5 -j+.5], 'k:')
                        end
                    end
                end
            end
        end

        for st1=1:length(R.St1)
            S1 = R.St1{st1}';
            SL1 = R.spikeLabel1{st1};
            if ~isempty(P.restrict2Time)
                idx = S1>P.restrict2Time(1) & S1<P.restrict2Time(2);
                S1 = S1(idx);
                SL1 = SL1(idx);
            end
            plot(S1, st1*ones(1, length(S1))-.5, '.', 'color', mysort.plot.vectorColor(st1), 'markersize', 14);
            for i=1:length(S1)
                if any(SL1(i) == P.mark)
                    plot(S1(i), st1-.5, 'rd', 'markersize', 10, 'linewidth', 2);
                end
                if ~any(SL1(i) == P.doNotLabel)
                    text(S1(i), st1-.2, defs.labelID2String{SL1(i)});
                end
            end
        end
        for st2=1:length(R.St2)
            S2 = R.St2{st2}';
            SL2 = R.spikeLabel2{st2};
            if ~isempty(P.restrict2Time)
                idx = S2>P.restrict2Time(1) & S2<P.restrict2Time(2);
                S2 = S2(idx);
                SL2 = SL2(idx);
            end            
            plot(S2, -st2*ones(1, length(S2))+.5, '.', 'color', mysort.plot.vectorColor(st2+length(R.St1)), 'markersize', 14);
            idx = ismember(SL2, P.mark);
            if ~isempty(idx) && any(idx)
                plot(S2(idx),repmat(-st2+.5, sum(idx),1), 'rd', 'markersize', 10, 'linewidth', 2);
            end
            idx = ~ismember(SL2, P.doNotLabel);
            if ~isempty(idx) && any(idx)
                labs = defs.labelID2String(SL2(idx));
                text(S2(idx),repmat(-st2+.2, sum(idx), 1), labs);
            end
        end    

        axis tight
        xlim = get(gca, 'xlim');
        plot(xlim +[-10 10], [0 0], 'k-', 'linewidth', 2);
        axis tight
        set(gca, 'ylim', [-length(R.St2) max(1, length(R.St1))]); % make sure this is not [0 0] !
    end
end