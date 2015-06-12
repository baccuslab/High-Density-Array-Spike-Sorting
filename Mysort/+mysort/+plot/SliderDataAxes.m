classdef SliderDataAxes < guiutil.SliderAxes
    properties
        % internals
        timeMode
        timeMultiplicators
        dataSources
        minmax
        plotSortingIndicesPerDataSource
        dataColors
        bPlotSortings
        channelSpacers
        maxTemplateChannelsToPlot
        templateChannelThreshold
        zeroLineStyle
        plottedGdfs
        
        controlledAxesHandles
    end
 
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = SliderDataAxes(dataSources, varargin)
            P.timeIn = 'samples'; % 'sec', 'ms'
            P.minWindowWidth = 100;
            P.maxWindowWidth = 10000;
            P.dataColors = {};
            P.plotSortings = 1;
            P.channelSpacers = [];
            P.channelSpacer = [];
            P.maxTemplateChannels = 25;
            P.templateChannelThreshold =5;
            P.zeroLineStyle = {};
            P.controlledAxesHandles = [];
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);               
            self = self@guiutil.SliderAxes(uP{:}, 'name', 'time', 'displayName', 'samples');
            if ~isempty(P.channelSpacer)
                P.channelSpacers = P.channelSpacer;
            end
            self.bPlotSortings = P.plotSortings;
            self.setMode(P.timeIn);
            self.minmax = [];
            self.windowWidthMax = P.maxWindowWidth;
            self.windowWidthMin = P.minWindowWidth;
            if ischar(P.dataColors)
                P.dataColors = {P.dataColors};
            end
            self.dataColors = P.dataColors;
            self.channelSpacers = P.channelSpacers;
            self.maxTemplateChannelsToPlot = P.maxTemplateChannels;
            self.templateChannelThreshold = P.templateChannelThreshold;
            self.zeroLineStyle = P.zeroLineStyle;
            self.makeLayout();
            self.setDataSource(dataSources);
            self.bAvoidReentrant = false;
            self.controlledAxesHandles = P.controlledAxesHandles;
            set(self.Axes, 'ButtonDownFcn', @self.CBDataClick);
        end
        %------------------------------------------------------------------
        function CBDataClick(self, a, b)
            point = get(gca, 'CurrentPoint');
            
            self.plotSortingIDs();
        end
        %------------------------------------------------------------------
        function plotSortingIDs(self)
            if isempty(self.plottedGdfs)
                return
            end
            for s=1:length(self.plottedGdfs)
                for i=1:size(self.plottedGdfs{s},1)
                    tM = self.timeMultiplicators(s);
                    id = self.plottedGdfs{s}(i,1);
                    time = self.plottedGdfs{s}(i,2)*tM;
                    text(time, -id, num2str(id),...
                        'parent', self.Axes(s), 'BackgroundColor', 'w');
                end
            end
        end
        %------------------------------------------------------------------
        function setMode(self, timeMode)
            self.timeMode = timeMode;
            if ~isempty(self.dataSources)
                self.setModeAndDataSource();
            end
        end
        %------------------------------------------------------------------
        function setDataSource(self, ds)
            ds = mysort.ds.checkDataSources(ds);
            if length(ds) ~= length(self.Axes)
                self.makeLayout(length(ds));
            end
            self.bAvoidReentrant = false;
            self.dataSources = ds;
            if ~mysort.ds.dataSourceIsEmpty(ds)
                for i=1:length(self.dataSources)
                    if isfield(self.dataSources{i}, 'activeSpikeSortingIdx') && ...
                       ~isempty(self.dataSources{i}.activeSpikeSortingIdx)
                        self.plotSortingIndicesPerDataSource{i} = {self.dataSources{i}.activeSpikeSortingIdx};
                    else
                        self.plotSortingIndicesPerDataSource{i} = {1};
                    end
                end
                if isempty(self.channelSpacers)
                    self.channelSpacers = zeros(1, length(self.dataSources));
                end
                self.setModeAndDataSources();
                self.update();
            else
                self.plotSortingIndicesPerDataSource = {};
                self.reset();
            end
        end
        %------------------------------------------------------------------
        function setModeAndDataSources(self)
            p = self.SliderPositionEditProperties.getProperty('time');
            w = self.WindowEditProperties.getProperty('window');
            self.minmax = [];
            if strcmp(self.timeMode, 'ms')
                p.setDisplayName('ms');
                p.set(1);                
                xlabel(self.Axes(end), 'time in ms');
                w.set(100);  
                self.windowWidthMin = 3;
                self.windowWidthMax = 1000;
            elseif strcmp(self.timeMode, 'sec') || strcmp(self.timeMode, 's')
                p.setDisplayName('sec');
                p.set(1);                
                xlabel(self.Axes(end), 'time in sec');
                w.set(.1); 
                self.windowWidthMin = .003;
                self.windowWidthMax = 10;
            elseif strcmp(self.timeMode, 'samples')
                p.setDisplayName('samples');
                p.set(1);                
                xlabel(self.Axes(end), 'time in samples');
                w.set(1000);
                self.windowWidthMin = 100;
                self.windowWidthMax = 100000;
            else 
                error('unknown time mode');
            end 
            for i=1:length(self.dataSources)
                ds = self.dataSources{i};
                sps = ds.getSamplesPerSecond();
                x = ds(1,1);
                assert(isfloat(x) || isinteger(x), 'Data source must provide float or int!')
                title(self.Axes(i), ds.name);
                if strcmp(self.timeMode, 'ms')
                    assert(~isempty(sps), 'If ms timeMode is selected all datasources must have a valid sample rate!');
                    self.timeMultiplicators(i) = 1000/sps;         
                elseif strcmp(self.timeMode, 'sec') || strcmp(self.timeMode, 's')
                    assert(~isempty(sps), 'If sec timeMode is selected all datasources must have a valid sample rate!');
                    self.timeMultiplicators(i) = 1/sps;        
                elseif strcmp(self.timeMode, 'samples')
                    self.timeMultiplicators(i) = 1;             
                else
                    error('unknown display type');
                end       
            end  
            if length(self.bPlotSortings) < length(self.dataSources)
                self.bPlotSortings = repmat(self.bPlotSortings, length(self.dataSources), 1);
            end
   
            % calculate min slider step with .5 overlap of windows
            L = size(self.dataSources{1},1)*self.timeMultiplicators(1);
            steps = min(1/(2*L/w.get()),1);
            % this will trigger the update!
            set(self.Slider, 'Min', 0, 'Max', L, 'Value', 0, 'SliderStep',[steps .05]);            
        end
        %------------------------------------------------------------------
        function update(self)
%             fprintf('Update?');
            if ~self.bIsVisible
%                 fprintf(' No.\n');
                return
            end            
%             fprintf(' Yes.\n');
            if isempty(self.dataSources)
                self.reset();
                return
            end
            self.bAvoidReentrant = true;
%             caller_name = mysort.util.getCallerName();
%             disp(['update - ' caller_name])
            w = self.WindowEditProperties.getProperty('window');
            w.set(min(max(w.get(), self.windowWidthMin), self.windowWidthMax));
            wW = w.get();
            for i=1:length(self.dataSources)
                ds = self.dataSources{i};
                tM = self.timeMultiplicators(i);
                chSp = self.channelSpacers(i);
                ct = max(1, round(self.sliderPosition/tM));
                t2 = min(ct + round(wW/tM), ds.size(1));
                t1 = max(1, t2-round(wW/tM));
                X = ds(t1:t2, :);
                %X = X-repmat(mean(X,2), 1, size(X,2));
                mima = [min(X(:)) max(X(:))]; 
                mima(2) = mima(2)+size(X,2)*chSp;
                if isempty(self.minmax) || size(self.minmax,1) < i
                    self.minmax(i,:) = mima;
                else
                    self.minmax(i,:) = [min(self.minmax(i,1), mima(1)) max(self.minmax(i,2), mima(2))];
                end
                cla(self.Axes(i))
                mysort.plot.zoomReset(self.Axes(i));
                set(self.Axes(i), 'nextplot', 'add');
                
                if ~isempty(self.zeroLineStyle) && ~isempty(self.zeroLineStyle{i})
                    plot(self.Axes(i), [t1 t2]*tM, [0 0], self.zeroLineStyle{i});
                end
                
                if length(self.dataColors) < i || isempty(self.dataColors{i})
                    for ch=1:size(X,2)
                        plot(self.Axes(i), (t1:t2)*tM, (ch-1)*chSp+X(:,ch), 'color', mysort.plot.vectorColor(ch));
                    end
                else
                    for ch=1:size(X,2)
                        plot(self.Axes(i), (t1:t2)*tM, (ch-1)*chSp+X(:,ch), 'color', self.dataColors{i});
                    end                    
                end

                if self.bPlotSortings(i) && ds.hasSpikeSorting()
                    S = ds.getActiveSpikeSorting();
%                     ME = ds.getMultiElectrode();
                    if ~isempty(S)
                        gdf = S.getGdf(t1, t2);
                        self.plottedGdfs{i} = gdf;
                        if ~isempty(gdf)
                            unitNames = unique(gdf(:,1));
                            if isa(S, 'mysort.spiketrain.SpikeSortingContainer')

                                T = S.getTemplateWaveforms(unitNames);
                                T = mysort.wf.pruneTemplates(T, ...
                                    'maxChannels', self.maxTemplateChannelsToPlot,...
                                    'minChannels', 1, ...
                                    'absThreshold', self.templateChannelThreshold,...
                                    'setInvalidChannelsTo', nan);
                                cutLeft = S.getTemplateCutLeft();
                                if S.isGroundTruth
                                    plot_mode = 'groundtruth';
                                else
                                    plot_mode = 'normal';
                                end
                                mysort.plot.templateSpikeTrain(T, gdf,...
                                    'timeMultiplicator', tM, ...
                                    'channelSpacer', chSp, ...
                                    'mode', plot_mode,...
                                    'AxesHandle', self.Axes(i),... 
                                    'cutLeft', cutLeft,...
                                    'T_gdf_idx2id', unitNames);%'spacer', 0,...
                                self.minmax(i,:) = [min(self.minmax(i,1), min(T(:))) max(self.minmax(i,2), max(T(:)))];
                            else
                                us = unique(gdf(:,1));
                                dY = diff(self.minmax(i,:));
                                offsetY = dY/300;
                                for uidx = 1:length(us)
                                    u = us(uidx);
                                    timepoints = gdf(gdf(:,1) == u, 2);
                                    %mysort.plot.markerForClass(u)
                                    plot(self.Axes(i), timepoints, offsetY*(u-1)*ones(1,length(timepoints))-1, ...
                                        'marker', 'x',...
                                        'color', mysort.plot.vectorColor(u),...
                                        'linewidth', 2, ...
                                        'LineStyle', 'none',...
                                        'markersize', 10);
                                end
                            end
                        end
                    end
                end
                set(self.Axes(i), 'Xlim', [t1 t2]*tM, ...
                                  'Ylim', self.minmax(i,:));      
            end
            if ~isempty(self.controlledAxesHandles)
                set(self.controlledAxesHandles, 'Xlim', [t1 t2]*tM);      
            end            
%             self.bAvoidReentrant = false;
        end
        %------------------------------------------------------------------
        function reset(self)
            self.dataSources = {};
            self.timeMultiplicators(1) = 1;
            p = self.SliderPositionEditProperties.getProperty('time');
            p.setDisplayName('samples');
            p.set(1);
            set(self.Slider, 'Min',0,'Max',1,'Value',.5,'SliderStep',[0 1]);
            cla(self.Axes);
            self.sliderPosition = 1;
            self.bAvoidReentrant = false;
        end        
    end
end