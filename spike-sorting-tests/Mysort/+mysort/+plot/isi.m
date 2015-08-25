
function RET = isi(spiketrains_or_gdf, varargin)
    P.figure = 1;
    P.figureHandle= -1;
    P.isHist = 0;
    P.binning = [];
    P.counts = [];
    P.IDs = [];
    P.title = '';
    P.figureTitle = [];
    P.printPercViolation = [];
    P.print2ms = 0;
    P.axesHandles = [];
    P.plotYLabels = 1;   
    P.edges = 0:32:20*32;
    P.srate = 32000;
    P.showIDs = 1;
    P.barParams = {'FaceColor', 'b'};
    P = mysort.util.parseInputs(P, 'isi', varargin);
    
    [spikeTrains P.IDs] = mysort.util.checkSpikeTrainsOrGDF(spiketrains_or_gdf, P);
    numHists = length(P.IDs);
    if isempty(numHists) || numHists == 0
        return
    end
    RET.fh = 0;
    if P.figureHandle > -1
        RET.fh = P.figureHandle;
    elseif P.figure
        RET.fh = mysort.plot.figure();
        mysort.plot.figureName('ISI');
    end  
    
    if ~P.isHist
        numTrains = numHists;
        spikeTrainsHists = zeros(numTrains, length(P.edges));
        counts = zeros(numTrains,1);
        for i=1:numTrains
            myTrain = spikeTrains{P.IDs(i)};
            if length(myTrain)>1
                spikeTrainsHists(i,:) = histc(diff(myTrain), P.edges);
            end
            counts(i) = length(myTrain);
        end
        if isempty(P.counts)
            P.counts = counts;
        end
    end
    spikeTrains = spikeTrainsHists;
    
    % downwards compatibility:
    if ~isempty(P.print2ms) && P.print2ms>0
        P.printPercViolation = 2;
    end
    
    if isempty(P.axesHandles)
        P.axesHandles = mysort.plot.subplots(numHists, 'figureTitle', P.figureTitle);
    end
    
    if isempty(P.binning)
        P.binning = (P.edges/P.srate)*1000;
    end
    if isempty(P.counts)
        P.counts = sum(spikeTrains, 2);
    end
    
    for i=1:numHists
        set(P.axesHandles(i), 'nextPlot', 'replace'); 
        bar(P.axesHandles(i), P.binning+(P.binning(2)-P.binning(1))/2, spikeTrains(i,:), P.barParams{:});
        set(P.axesHandles(i), 'nextPlot', 'add'); 
        
        str = '';
        if ~isempty(P.IDs) && P.showIDs
            str = [str 'ID: ' num2str(P.IDs(i)) ' '];
            plot(P.axesHandles(i), P.binning(end-1), max(spikeTrains(i,:)/2), '.', ...
                'color', mysort.plot.vectorColor(P.IDs(i)), 'markerSize', 25);                   
        end
        
        str = [str '#: ' num2str(P.counts(i)) ' '];
        
        if ~isempty(P.printPercViolation)
            p2ms(i) = 100*sum(spikeTrains(i,1:find(P.binning>P.printPercViolation,1)-1))/P.counts(i);
            str = [str num2str(P.printPercViolation) 'ms: ' num2str(round(p2ms(i)))];
        end
        if P.plotYLabels
            ylabel(P.axesHandles(i), str);
        end
        if i==1 && ~isempty(P.title)
            title(P.axesHandles(i), P.title);
        end
%         if i~=numHists
%             set(gca, 'xTickLabel',[]);
%         end
        if i==numHists
            xlabel(P.axesHandles(i), 'Time in ms')
        end
        title(P.axesHandles(i), sprintf('Class %d', P.IDs(i)))
        set(P.axesHandles(i), 'xlim', [0 P.binning(end)]);
    end        
    if nargout >0
        RET.p3ms = p3ms;
    end
end