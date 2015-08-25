
function sorting(X, gdf, T, varargin)
    P.figureHandle = [];
    P.axesHandle = [];
    P.sampleOffset = 0;
    P.srate = [];
    P.Y = [];
    P.spacer = [];
    P.cutleft = 0;
    P.gtGdf = [];
    P.template_IDs = [];
    P.restrictToChannels = [];
    P = mysort.util.parseInputs(P, 'sorting', varargin);
    
    if isempty(P.axesHandle)
        if isempty(P.figureHandle)
            P.figureHandle = mysort.plot.figure('color','w');
            mysort.plot.figureName('Sorting');
        end
        P.axesHandle = axes();
    end
%     axes(P.axesHandle);
    if isempty(P.restrictToChannels)
        P.restrictToChannels = 1:size(X,1);
    end
    P.spacer = mysort.plot.mc(X(P.restrictToChannels,:), 'figure', 0, 'color', {'k'},...
                'spacer', P.spacer, 'sampleOffset', P.sampleOffset,...
                'axesHandle', P.axesHandle);
    set(P.axesHandle, 'nextplot', 'add'); 

    if isempty(P.template_IDs)
        P.template_IDs = 1:size(T,3);
    end
    
    T = T(:,P.restrictToChannels,:);
    
    if ~isempty(gdf)
        gdf(:,2) = gdf(:,2) - P.cutleft;
        mysort.plot.templateSpikeTrain(T, gdf, 'channelSpacer', P.spacer, 'mode', 'normal', 'AxesHandle', P.axesHandle); %, 'template_IDs', P.template_IDs)
    end    
    if ~isempty(P.gtGdf)
        mysort.plot.templateSpikeTrain(T, P.gtGdf, 'channelSpacer', P.spacer, 'mode', 'gt', 'AxesHandle', P.axesHandle); %, 'template_IDs', P.template_IDs)
    end
    axis tight
    