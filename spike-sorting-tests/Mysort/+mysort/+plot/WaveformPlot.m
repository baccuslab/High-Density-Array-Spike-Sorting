classdef WaveformPlot < guiutil.SingleAxes
    properties
        MultiElectrode
        wfs
        wfs_ids
        
        % gui elements
        PropertyList
        GuiProperties
        
        % gui internals
        disableCallback
        
        m_nPlotMode
        m_nMaxPlotWaveforms
        m_nMaxPlotChannels
        
        plotControls
        markedElectrodeNumbers
        electrodeMarkerHandle
        waveformHandle
        templateHandle
        
        % externals
        waveformColor
        templateColor
        
    end
    
 
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = WaveformPlot(wfs, ME, varargin)
            %P.plotMode = 'horizontal'; % 'vertical', '2d'
            P.template = [];
            P.waveformColor = [.5 .5 .5];
            P.templateColor = [.1 .1 .9];
            P.unit = [];
            P.waveformIDs = [];
            P.plotTemplate = false;
            P.plotControls = false;
            P.maxPlotWaveforms = 1000;
            P.maxPlotChannels = 10;    
            P.plotMode = 'horizontal'; % 'vertical', '2d'            
           
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);       
            self = self@guiutil.SingleAxes(uP{:});
            if ~isempty(P.unit)
                P.templateColor = mysort.plot.vectorColor(P.unit);
            end
            self.wfs_ids = P.waveformIDs;
            self.m_nMaxPlotWaveforms = P.maxPlotWaveforms;
            self.m_nMaxPlotChannels = P.maxPlotChannels;
            self.plotControls = P.plotControls;
            self.waveformColor = P.waveformColor;
            self.templateColor = P.templateColor;
            self.m_nPlotMode = P.plotMode;        
            self.makeLayout();
            self.setWfAndMultiElectrode(wfs, ME);
            
            self.disableCallback = 1;
            self.GuiProperties.getProperty('plotMode').set(self.m_nPlotMode);
            self.disableCallback = 0;            
        end
        %------------------------------------------------------------------
        function setWfAndMultiElectrode(self, wfs, ME)
            self.reset();
            if isnumeric(ME)
                ME = mysort.ds.MultiElectrode(ME, 1:size(ME,1));
            end
            self.MultiElectrode = ME;
               
            self.wfs = wfs;
            if ~isempty(wfs)
                self.update();
            end
        end
        %------------------------------------------------------------------
        function plot(self, wfs, ME)
            setWfAndMultiElectrode(self, wfs, ME);
        end
        %------------------------------------------------------------------
        function reset(self)
            cla(self.Axes);
            mysort.plot.zoomReset(self.Axes);
        end
        %------------------------------------------------------------------
        function makeLayout(self)
            p = self.Parent;
            
            self.initPropertyList();
            
            if self.plotControls
                mainLayout = uiextras.HBox( 'Parent', p, 'Spacing', 3);
                    wfbox = uiextras.BoxPanel( ...
                       'Parent', mainLayout, ...
                       'Title', 'Waveforms');
                        if isempty(self.Axes)
                            self.Axes = axes('Parent', wfbox);
                        end                    
                    controlbox = uiextras.BoxPanel( ...
                       'Parent', mainLayout, ...
                       'Title', 'Controls');
                        % create property interface
                        self.GuiProperties = guiutil.GuiProperties(controlbox, self.PropertyList,...
                                'callback', [],...
                                'loadSaveButton', 0, ...
                                'labelArrangement', 'vertical');  

                set(mainLayout, 'Sizes', [-1, 100]);
            else
                self.GuiProperties = guiutil.PropertyList(self.PropertyList);
                mainLayout = uiextras.HBox( 'Parent', p, 'Spacing', 3);
                    if isempty(self.Axes)
                        self.Axes = axes('Parent', mainLayout);
                    end                    
                set(mainLayout, 'Sizes', [-1]);                
            end
        end
        %------------------------------------------------------------------
        function update(self)
            self.reset();
            if isempty(self.MultiElectrode) || isempty(self.wfs)
                return
            end            
            self.m_nPlotMode = self.GuiProperties.getProperty('plotMode').get();
            self.plotWaveforms();
            self.plotElectrodeMarker();
        end

        %------------------------------------------------------------------
        function plotWaveforms(self, wfs, EP, varargin)
            P.templateColor = self.templateColor;
            P.waveformColor = self.waveformColor;
            P.plotTemplate = self.GuiProperties.getProperty('plotTemplate').get();
            P.template = [];
            P.templatePlotArgs = {'-', 'linewidth', 2};
            P.plotArgs = {'-', 'linewidth', 1};
            P = mysort.util.parseInputs(P, varargin, 'error');
            if nargin < 3
                EP = self.MultiElectrode.electrodePositions;
            end
            if nargin < 2
                wfs = self.wfs;
            end
            self.waveformHandle = mysort.plot.waveforms2D(wfs, EP,...
                'AxesHandle', self.Axes, ...
                'plotArgs', [P.plotArgs {'color', P.waveformColor}]);
            if P.plotTemplate 
                self.plotTemplate(P.template, wfs, EP, P.templatePlotArgs{:}, 'color', P.templateColor);
            end
        end
        %------------------------------------------------------------------
        function plotTemplate(self, template, wfs, EP, varargin)  
            if isempty(template)
                nC = size(wfs,2);
                template = mysort.wf.v2t(median(mysort.wf.t2v(wfs),1), nC);
            end
            self.templateHandle =  mysort.plot.waveforms2D(template, EP,...
                'AxesHandle', self.Axes, ...
                'plotArgs', varargin);
        end
        %------------------------------------------------------------------
        function plotElectrodeMarker(self)
            self.removeElectrodeMarker();
            idx = self.MultiElectrode.getElIdx4ElNumber(self.markedElectrodeNumbers);
            if isempty(idx)
                return
            end
            EP = self.MultiElectrode.electrodePositions(idx,:);
            EN = self.MultiElectrode.electrodeNumbers(idx);
            self.electrodeMarkerHandle = plot(self.Axes, EP(:,1), EP(:,2), '*', 'color', 'm', 'markersize', 8, 'linewidth', 2);
            for i=1:length(idx)
                self.electrodeMarkerHandle(i+1) = text(EP(i,1), EP(i,2), num2str(EN(i)), 'Parent', self.Axes);
            end
            set(self.electrodeMarkerHandle, 'HitTest', 'off');
        end
        %------------------------------------------------------------------
        function removeElectrodeMarker(self)
            if ishandle(self.electrodeMarkerHandle)
                try
                    delete(self.electrodeMarkerHandle)
                catch
                end
            end
        end
        %------------------------------------------------------------------
        function setMarkedElectrodeNumbers(self, elNumbers)
            self.markedElectrodeNumbers = elNumbers;
            self.plotElectrodeMarker();
        end   
        %------------------------------------------------------------------
        function CBPropertyChange(self, a)
            if self.disableCallback==1
                return
            end
            self.update();
        end
        %------------------------------------------------------------------
        function initPropertyList(self)
            str = {};
            if ~isempty(self.wfs_ids)
                str = arrayfun(@(x) num2str(x), unique(self.wfs_ids), 'UniformOutput', false);
            end
            props = {'Unit', {'all' str{:}}
                     'plotMode', {'horizontal', 'vertical', '2D'}
                     'plotMaxChan', self.m_nMaxPlotChannels
                     'plotNWfs', self.m_nMaxPlotWaveforms
                     'plotRandWfs', false
                     'plotElNumbers', true
                     'plotTemplate', true};            
            self.PropertyList = guiutil.Property.empty();
            cb = @self.CBPropertyChange;
            for i=1:size(props,1)
                self.PropertyList(end+1) = guiutil.Property(props{i,:}, 'changeCallback', cb);
            end       
        end
    end
end