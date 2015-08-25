
function [varargout] = mc(varargin)
    import mysort.*
    % X - Multichannel Data, first Argument!
    % 'start'  - startpoint in timesteps
    % 'end'    - endpoint in timesteps
    % 'events' - matrix of events (1. column names, 2. timestep)
    % 'srate'  - samplerate

    [X, P, D] = parseInputs(varargin);
    
    assert(~isempty(P.sampleOffset), 'sampleOffset must not be empty!');
    
    if ~isempty(P.axesHandle)
        axes(P.axesHandle);
        if isempty(P.figureHandle)
            P.figureHandle = gcf; 
        end
    elseif ~isempty(P.figureHandle)
        mysort.plot.figure(P.figureHandle);
        clf
        mysort.plot.figureName('mc');
    elseif(P.figure==1)
        P.figureHandle = mysort.plot.figure;
        mysort.plot.figureName('mc');
        a = axes;
    else
        P.figureHandle = gcf;
    end
    
    
    
    if isempty(X)
        varargout{1} = 0;
        varargout{2} = P.figureHandle;
%        warning('X is empty in PlotMC!');
        return
    end

    if P.normalizeData == 1
        X = X-repmat(mean(X,2), 1, size(X,2));
        X = X./repmat(std(X,0,2), 1, size(X,2));
    end

    hold on
    if isempty(P.spacer) || (P.spacer == -1)
        if D.channelnum==1
            M = max(X) + abs(min(X));
        else
            M = max(max(X(1:(end-1),:)))+abs(min(min(X(2:end,:))));
        end
    else
        M = P.spacer;
    end
    %MA = max(X');
    %MI = [abs(min(X')) 0];
    offset = 0;
    P.sampleOffset = P.sampleOffset*P.srate;
    if size(P.xticks) == [0 0]
        xticks = [1:size(X(1,:),2)]*P.srate;
    else
        xticks = P.xticks;
    end
    
    if size(xticks,1) == 1 && D.channelnum >1
        xticks = repmat(xticks, D.channelnum,1);
    end
    
    if size(P.colors) == [1 1];
       P.colors = repmat(P.colors, D.channelnum, 1);
    end
    for c = 1:D.channelnum
        if P.eventsfirst==0
            %plot(xticks, X(c,:)+offset,P.color);
            plot(P.sampleOffset+ xticks(c,:), P.yOffset+X(c,:)+(c-1)*M,P.colors{c},'LineSize',P.linewidth);
            hold on
            if P.showZeroLine
                plot(P.sampleOffset+ [xticks(c,1) xticks(c,end)], P.yOffset+[0 0], ':k');
            end
        end
        %offset = offset + MA(c)+MI(c+1);
    end
    if D.channelnum == 1
        offset = 0;
    end

    if ~isempty(P.events)
        factor = 1;
        if ~P.eventsInTime 
            factor = P.srate;
        end        
        if iscell(P.events)
            P.events = cellEvents2StructEvents(P.events);
        end
        if isnumeric(P.events)
            P.events = matEvents2StructEvents(P.events);
        end
        ymi = min(X(1,:));
        yma = 1.1*max(X(end,:));
        for i=1:length(P.events)

            xx = [P.events(i).sample P.events(i).sample]*factor;
            yy = P.yOffset+[ymi M*(D.channelnum-1)+yma];
            plot(xx,yy,'Color', mysort.plot.vectorColor(P.events(i).channel),'LineStyle',':','Linewidth',2);
            text(xx(1),yy(2),P.events(i).name);
        end
    end

    if ~isempty(P.eventsTime)
        for e = 1:size(P.eventsTime,1)
            %l=line([P.eventsTime(e,end) P.eventsTime(e,end)], [min(X(1,:)) offset + MA(end)]);
            l=line(P.sampleOffset+[P.eventsTime(e,end) P.eventsTime(e,end)], P.yOffset+[min(X(1,:)) M*(D.channelnum-1)+1.1*max(X(end,:))]);
            set(l, 'Color', mysort.plot.vectorColor(3));
        end
    end

    if P.eventsfirst==1
        %offset = 0;
        for c = 1:D.channelnum
            %plot(xticks, X(c,:)+offset,P.color);
    %        try
            if ischar(P.colors{c})
                plot(P.sampleOffset+ xticks(c,:), P.yOffset+X(c,:)+(c-1)*M,...
                    P.colors{c},'LineWidth',P.linewidth,'LineStyle',P.lineStyle);
            else
                plot(P.sampleOffset+ xticks(c,:), P.yOffset+X(c,:)+(c-1)*M,...
                    'color', P.colors{c},'LineWidth',P.linewidth,'LineStyle',P.lineStyle);
            end
            hold on
            if P.showZeroLine
                plot(P.sampleOffset+ [xticks(c,1) xticks(c,end)], P.yOffset+(c-1)*M+[0 0], ':k');           
            end
    %        catch
    %            display('Error in PlotMC');
    %        end
            %offset = offset + MA(c)+MI(c+1);
        end
    end

    if ~P.yLabels
        set(gca, 'YTickLabel',[]);
    end



    if isfield(P,'threshold') && P.threshold < inf;
        for i=1:length(P.threshold)
            thresh = P.threshold(length(P.threshold)-i+1);
            if ~isnumeric(P.thresholdColor)
                plot(P.sampleOffset+ [xticks(i,1) xticks(i,end)], P.yOffset+(i-1)*M+[thresh thresh], P.thresholdColor);
            else
                plot(P.sampleOffset+ [xticks(i,1) xticks(i,end)], P.yOffset+(i-1)*M+[thresh thresh], 'color', P.thresholdColor);
            end
        end
    end

    if P.marker ~= 0
        plot(P.sampleOffset+ min(xticks(:,1)), P.yOffset+P.markerHeight, '.g', 'MarkerSize',16);
    end
    
    if ~isempty(P.epochs)
        for i=1:size(P.epochs,1)
            plot(P.epochs(i,:), zeros(1,2), '-r', 'lineWidth', 3);
        end
    end

    % Set Data Cursor precision not to use the stupid standard format
    mysort.plot.dataCursorMode(P.figureHandle);
    
    if nargout > 0
        varargout{1} = M;
        varargout{2} = P.figureHandle;
    end    



    % -------------------------------------------------------------------------
    function E = cellEvents2StructEvents(CE)
        count = 0;
        for cIdx=1:length(CE)
            for ii=1:length(CE{cIdx})
                count = count +1;
                E(count).sample  = CE{cIdx}(ii);
                E(count).channel = cIdx;
                E(count).name = num2str(cIdx);
            end
        end
    end
    % -------------------------------------------------------------------------
    function E = matEvents2StructEvents(ME)
        for mIdx=1:length(ME)
            E(mIdx).sample  = ME(mIdx,2);
            E(mIdx).channel = ME(mIdx,1);
            E(mIdx).name = num2str(ME(mIdx,1));
        end
    end
    % -------------------------------------------------------------------------
    function  [X, P, D]=parseInputs(in)
        X = squeeze(in{1});
        [nC,L] = size(X);
        if nC > 1000 && nC > L
            X = X';
            [nC,L] = size(X);
        end
        P.t1 = 1;
        P.t2 = size(X,2); %min(100000, size(X,2));
        P.srate = 1;
        P.epochs = [];
        P.events = [];
        P.eventsTime = [];
        P.eventsInTime = 0;
        P.eventColors = [0 1 0
                         0 1 1
                         1 0 0
                         1 1 0
                         0 0 1
                         0 0 0
                         1 0 1
                         .5 .5 .5
                         0  .5 .5
                         .5 0  .5
                         .5 .5 0
                         0  0  .5
                         0  .5 0
                         .5 0  0];
        P.eventTags = {};             
        P.colors = {'b'};
        P.color = P.colors;
        P.fig = 1;
        P.eventsfirst = 1;
        P.threshold = inf;
        P.thresholdColor = 'g';
        P.spacer = -1;
        P.xticks = [];
        P.linewidth = 1;
        P.lineStyle = '-';
        P.sampleOffset = 0;
        P.normalizeData = 0;
        P.yOffset = 0;
        P.yLabels = 1;
        P.figureHandle = [];
        P.marker = 0;
        P.markerHeight = 0;
        P.plotZeroLine = 0;
        P.figure = 1;
        P.axesHandle = [];
        P.showZeroLine = 0;
            for ii=2:2:size(in,2)
                if isfield(P, in{ii})
                    P.(in{ii}) = in{ii+1};
                else
                    warning(sprintf('Field not found in PlotMC: %s', in{ii}));
                end
            end
        P.colors = P.color;

        if(P.t2 == -1)
            P.t2 = L;
        end
        D.channelnum = nC;
        D.length = L;
        X = X([nC:-1:1],P.t1:P.t2);

    end
end
