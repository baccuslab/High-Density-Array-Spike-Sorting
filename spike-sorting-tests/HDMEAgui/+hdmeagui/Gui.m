classdef Gui < handle
    properties
        % gui elements
        Window
        mainLayout
            ControlTabPanel
                MEAElement
                Properties
            SortingControlLayout
                sortingControlButtons
                sortingTemplates
            ViewTabPanel
                DataInspector
                SpikeInspector
                SortingInspector
                TemplateMatchingInspector
                WaveformInspector
        
        % internals
        TabGuiElems
        
        selectedSortingWaveformManager
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function gui = Gui()
            gui.Window = mysort.plot.figure('w', 1200, 'h', 800, 'name', 'HDMEA GUI');
            gui.mainLayout = uiextras.HBoxFlex( 'Parent', gui.Window, 'Spacing', 3 );

                % Create the panels
                gui.ControlTabPanel = uiextras.TabPanel( 'Parent', gui.mainLayout,...
                    'Padding', 5 );
                gui.ViewTabPanel = uiextras.TabPanel( 'Parent', gui.mainLayout,...
                    'Padding', 5, 'Callback', @gui.CBViewTabSelection, 'TabSize', 60);
                    gui.DataInspector = mysort.plot.SliderDataAxes([], ...
                        'Parent', gui.ViewTabPanel, 'visible', true, ...
                        'plotSortings', true, 'timeIn', 'sec',...
                        'channelSpacers', 50);
                    gui.SpikeInspector = hdmeagui.SpikeGui('Parent', gui.ViewTabPanel, 'visible', false);
                    gui.SortingInspector = hdmeagui.SortingGui('Parent', gui.ViewTabPanel,...
                        'visible', false, 'newSortingCallback', @gui.CBNewSorting);
                    gui.TemplateMatchingInspector = hdmeagui.TemplateMatchingGui('Parent', gui.ViewTabPanel,...
                        'visible', false, 'newSortingCallback', @gui.CBNewSorting);
                    gui.WaveformInspector = mysort.plot.MultiSessionWaveformManagerPlot([], ...
                                'maxPlotWaveforms', 100, ...
                                'maxPlotChannels', 5,...
                                'plotMode', '2D',...
                                'plotControls', 1,...
                                'Parent', gui.ViewTabPanel);                  
                gui.ViewTabPanel.TabNames = {'Data', 'Spikes', 'Sorting', 'TemplateM', 'Wfs'};
                gui.ViewTabPanel.SelectedChild = 1; 
            % Adjust the main layout
            set( gui.mainLayout, 'Sizes', [-1,-5]);
            gui.TabGuiElems = {gui.DataInspector, gui.SpikeInspector, gui.SortingInspector,...
                gui.TemplateMatchingInspector, gui.WaveformInspector};
            % fill the control tab
            gui.makeControlTab();
        end
        
        %------------------------------------------------------------------
        function makeControlTab(gui)
            gui.MEAElement = mysort.mea.H5FileGuiElement(gui.ControlTabPanel,...
                'hasOpenButton', 1,...
                'openCallback', @gui.CBMEAOpenButton,...
                'selectionCallback', @gui.CBMEASessionSelection,...
                'sortingSelectionCallback', @gui.CBMEASortingSelection,...
                'templateSelectionCallback', @gui.CBMEATemplateSelection);
            
            % load property values
            if exist('props.mat', 'file')
                D = load('props.mat');
                props = D.props;
            else
                props = {'cutTf', 35
                         'cutLeft', 10
                         'showTf', 35
                         'showLeft', 10
                         'testBool', true};
            end
            cb = @gui.CBPropertyChange;
           
            % create property interface
            gui.Properties = guiutil.GuiProperties(gui.ControlTabPanel, props,...
                'callback', cb,...
                'loadSaveButton', 1);
            
            gui.ControlTabPanel.TabNames = {'MEA File', 'Gui Props'};
            gui.ControlTabPanel.SelectedChild = 1; 
        end

        
        %------------------------------------------------------------------
        function CBMEAOpenButton(self)
            disp('file open')
            self.DataInspector.setDataSource(self.MEAElement.MEA);
            self.SpikeInspector.setDataSource(self.MEAElement.MEA);
            self.SortingInspector.setDataSource(self.MEAElement.MEA);
            self.TemplateMatchingInspector.setDataSource(self.MEAElement.MEA);
            self.selectedSortingWaveformManager = [];
        end
        %------------------------------------------------------------------
        function CBMEASessionSelection(self)
            self.DataInspector.setDataSource(self.MEAElement.MEA);
            self.SpikeInspector.setDataSource(self.MEAElement.MEA);
            self.SortingInspector.setDataSource(self.MEAElement.MEA);
            self.TemplateMatchingInspector.setDataSource(self.MEAElement.MEA);
            self.selectedSortingWaveformManager = [];
        end
        %------------------------------------------------------------------
        function CBMEASortingSelection(self, v)
            self.DataInspector.CBSetSliderPositionButton();
            self.SortingInspector.setSelectedSortingIndex(v);
            self.TemplateMatchingInspector.setSelectedSortingIndex(v);
            self.selectedSortingWaveformManager = [];
        end        
        %------------------------------------------------------------------
        function CBMEATemplateSelection(self, v)
            if isempty(self.selectedSortingWaveformManager)
                SC = self.MEAElement.SortingElement.selectedSortingContainer;
                ME = self.MEAElement.MEA.getAllSessionsMergedMultiElectrode();
                sidx = self.MEAElement.MEA.getSelectedSessionIdx();

                mGdf = SC.getMultiSessionGdf(sidx);
                nSp = size(mGdf,1);
                if nSp>0
                    self.selectedSortingWaveformManager = ...
                        mysort.wf.MultiSessionBufferedWfManager(self.MEAElement.MEA,...
                            mGdf(:,2), mGdf(:,1), mGdf(:,3), 20, 80, 1:nSp, ME);    
                end
            end
%             SC = self.MEAElement.SortingElement.selectedSortingContainer;
%             ME = self.MEAElement.MEA.getAllSessionsMergedMultiElectrode();
%             wfs = SC.TemplateManager.getWaveforms4MultiElectrode(ME);
%             wfmp = mysort.plot.WaveformPlot(wfs(:,:,1), ME, ...
%                 'plotTemplate', false, ...
%                 'plotControls', 0, ...
%                 'maxPlotChannels', 10, ...
%                 'waveformColor', mysort.plot.vectorColor(1));
%             set(wfmp.Axes, 'NextPlot', 'add');
%             for i=2:size(wfs,3)
%                 wfmp.plotWaveforms(wfs(:,:,i), ME.electrodePositions, ...
%                     'plotTemplate', false, 'waveformColor', mysort.plot.vectorColor(i));
%             end
            if isempty(v)
                self.WaveformInspector.setWfManager([]);
            else
                for i=1:length(v)
                    T = self.MEAElement.SortingElement.selectedSortingContainer.TemplateManager.TemplateList(v(i));
                    self.WaveformInspector.m_nRestrict2DetectionChannels = T.name;
                    c = mysort.plot.vectorColor(T.name);
                    self.WaveformInspector.templateColor = c;
                    self.WaveformInspector.setWfManager(self.selectedSortingWaveformManager);                
                    self.WaveformInspector.plotWaveforms([], T.MultiElectrode.electrodePositions, 'template', T.waveforms, 'templatePlotArgs', {'--'});
                end
            end
        end            
        %------------------------------------------------------------------
        function CBViewTabSelection(self, obj, ev)
            self.TabGuiElems{ev.SelectedChild}.setVisibility(true);
            self.TabGuiElems{ev.PreviousChild}.setVisibility(false);
        end
        %------------------------------------------------------------------
        function CBPropertyChange(self, p)
            disp(['property changed: ' p.getName()])
            disp(p.get())
        end        
        %------------------------------------------------------------------
        function CBNewSorting(self)
            self.MEAElement.update();
        end
    end
end

