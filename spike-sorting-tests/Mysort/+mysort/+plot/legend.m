function [h, lh] = legend(ax, plotArgs, plotNames, location)
    % produces a legend in axeshandle ax with the linetypes, colors and
    % names specified in varargin
    %
    % mysort.plot.legend(ah, {{'.r', 'linewidth', 2}, {'xb', 'linewidth', 2}, {'.k', 'linewidth', 2}}, {'T1', 'T2', 'T3'}, 'SouthEast')
    if ~ishandle(ax)
        if nargin == 3
            location = plotNames;
        else
            location = 'NorthEast';
        end
        plotNames = plotArgs;
        plotArgs = ax;
        ax = gca;
    elseif nargin == 3
        location = 'NorthEast';
    end
    
    % set all currently plotted elements off
    ch = get(ax, 'Children');
    for i=1:length(ch)
        set(get(get(ch(i),'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
    end

    % plot invisible objects representing the legend entries
    hold on
    for i=1:length(plotArgs)
        args = plotArgs{i};
        h = plot(ax,1,1,args{:});
        set(h, 'Visible', 'off', 'DisplayName', plotNames{i})
    end

    set(get(get(h,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');

    lh = legend(ax,'show');
    set(lh, 'location', location, 'box', 'off')