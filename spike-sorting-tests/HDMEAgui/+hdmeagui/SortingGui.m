classdef SortingGui < guiutil.GuiElement
    properties
        % gui elements
        FeaturePlot
        SpikeWaveformManagerPlot
        ChipWaveformPlot
        
        % gui internals
        AcceptTemplateButton
        NewSortingButton
        
        % internals
        dataSource
        preloadedData
        Features
        WaveformManager
        SelectedWaveformManager
        
        SelectedSortingIdx
        SpikeSortingContainer
        selectedTemplateIdx
        
        
        % callbacks 
        CBNewSorting
    end
 
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = SortingGui(varargin)
            P.dataSource = [];
            P.newSortingCallback = [];
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);
            self = self@guiutil.GuiElement(uP{:});
            self.makeLayout();
            self.dataSource = P.dataSource;
            self.CBNewSorting = P.newSortingCallback;
        end
        %------------------------------------------------------------------
        function makeLayout(self)
            p = self.Parent;            
            mainLayout = uiextras.HBox( 'Parent', p, 'Padding', 0 ); 
                viewLayout = uiextras.HBoxFlex( 'Parent', mainLayout,'Padding', 2 );
                    detailLayout = uiextras.VBoxFlex( 'Parent', viewLayout,'Padding', 2 );
                        controlButtonBox = uiextras.BoxPanel( ...
                               'Parent', detailLayout, ...
                               'Title', 'Controls');
                           controlButtonLayout = uiextras.HBox( 'Parent', controlButtonBox, 'Padding', 0 ); 
                               self.NewSortingButton = uicontrol( 'Style', 'PushButton', ...
                                       'Parent', controlButtonLayout, ...
                                       'String', 'NewSorting', ...
                                       'Callback', @self.CBNewSortingButton); 
                               self.AcceptTemplateButton = uicontrol( 'Style', 'PushButton', ...
                                       'Parent', controlButtonLayout, ...
                                       'String', 'Accept', ...
                                       'Callback', @self.CBAcceptTemplateButton); 
                               uiextras.Empty( 'Parent', controlButtonLayout )
                           set(controlButtonLayout, 'Sizes', [100, 100, -1]);
                        box = uiextras.BoxPanel( ...
                               'Parent', detailLayout, ...
                               'Title', 'Spike Selection');
                            someLayout = uiextras.VBoxFlex( 'Parent', box,'Padding', 2 );   
                                self.FeaturePlot = hdmeagui.FeaturePlot(hdmeagui.FeatureContainer([]), @self.CBfeatureSelection, 'Parent', someLayout);
                                self.SpikeWaveformManagerPlot = mysort.plot.MultiSessionWaveformManagerPlot([],...
                                    'Parent', someLayout, ...
                                    'plotMedianColor', [.9 .7 .0],...
                                    'maxPlotWaveforms', 200);
                            set(someLayout, 'Sizes', [-1,-1]);
                        selectionDetailsBox = uiextras.HBox( 'Parent', detailLayout,...
                        'Padding', 0 );
                            Box2 = uiextras.BoxPanel( ...
                               'Parent', selectionDetailsBox, ...
                               'Title', 'Box2');
                            Box3 = uiextras.BoxPanel( ...
                               'Parent', selectionDetailsBox, ...
                               'Title', 'Box3');
                        set( selectionDetailsBox, 'Sizes', [-2,-1]);
                    set( detailLayout, 'Sizes', [50, -10,-1]);
                    box = uiextras.BoxPanel( ...
                               'Parent', viewLayout, ...
                               'Title', 'Footprint');
                    self.ChipWaveformPlot = mysort.plot.WaveformPlot([], [],...
                        'maxPlotWaveforms', 1000, ...
                        'maxPlotChannels', 10000, 'Parent', box);
                set( viewLayout, 'Sizes', [-3,-1]);
             set(mainLayout, 'Sizes', [-1]);
        end
        %------------------------------------------------------------------
        function update(self)
            if ~self.bIsVisible || isempty(self.dataSource)
                return
            end
            self.preloadData();
            self.removeSpikesFromFeaturesThatAreTooCloseToAlreadySortedSpikes();
            self.FeaturePlot.setFeatureContainer(self.Features);      
        end  
        %------------------------------------------------------------------
        function self = removeSpikesFromFeaturesThatAreTooCloseToAlreadySortedSpikes(self)
            if isempty(self.SpikeSortingContainer)
                return
            end
            t1 = tic;  fprintf('Removing single electrode detected spikes that are explained by current sorting...');
            for i=1:length(self.Features.featureList)
                self.Features.featureList(i, 1).removeCloseEvents(self.SpikeSortingContainer.singleSessionGdfList, 5);
            end        
            t2 = toc(t1); fprintf(' done. (%.3f)\n', t2);
        end
        %------------------------------------------------------------------
        function self = setSelectedSortingIndex(self, v)
            if ~isempty(self.dataSource) && (isempty(self.SelectedSortingIdx) || v ~= self.SelectedSortingIdx)
                self.reset();
                self.SelectedSortingIdx = v;
                if ~isempty(self.dataSource.SpikeSortingContainers)
                    self.SpikeSortingContainer = self.dataSource.SpikeSortingContainers{v};
                end                
                self.update();                    
            end
        end        
        %------------------------------------------------------------------
        function CBNewSortingButton(self, a, b)
            prompt = {'Enter name for new spike sorting:'};
            dlg_title = 'Input for spike sorting name';
            num_lines = 1;
            def = {'no_name'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            if ~isempty(answer)     
                self.SpikeSortingContainer = self.dataSource.createNewSpikeSortingContainer(answer{1});
                if ~isempty(self.CBNewSorting)
                    self.CBNewSorting();
                end
            end                        
        end
        %------------------------------------------------------------------
        function checkSpikeSortingContainer(self)
            if isempty(self.dataSource)
                return
            end
            if ~isempty(self.dataSource.SpikeSortingContainers)
                return
            end
            self.CBNewSortingButton();
        end
        %------------------------------------------------------------------
        function reset(self)
            self.SpikeWaveformManagerPlot.reset();
            self.preloadedData = [];   
            self.SelectedSortingIdx = [];
            self.SpikeSortingContainer = [];
            self.selectedTemplateIdx = [];
        end  
        %------------------------------------------------------------------
        function CBfeatureSelection(self, F, thr)
            ids = F.getIDs4Threshold(thr);
            if ~isempty(ids)
                self.SelectedWaveformManager = self.WaveformManager.getWfManger4SubIDs(ids);
%                 disp(self.SelectedWaveformManager)
                self.SpikeWaveformManagerPlot.setWfManager(self.SelectedWaveformManager);
                self.ChipWaveformPlot.plot(self.SpikeWaveformManagerPlot.T', self.SpikeWaveformManagerPlot.ME);
                self.ChipWaveformPlot.setMarkedElectrodeNumbers(self.SpikeWaveformManagerPlot.plottedChannelNumbers);
                
                self.FeaturePlot.setMarkedFeatureElectrodes(self.SpikeWaveformManagerPlot.plottedChannelNumbers);
            else
                self.SpikeWaveformManagerPlot.reset();
            end
        end
        %------------------------------------------------------------------
        function preloadData(self)
            if ~isempty(self.preloadedData)
                return
            end
            t1 = tic;  fprintf('Loading single electrode detected Spikes...');
            D = self.dataSource.getSelectedSessionsDetectedSpikes();
            ME = self.dataSource.getSelectedSessionsMergedMultiElectrode();
            en = ME.getElectrodeNumbers();
            nC = ME.getNElectrodes();
%             ridx = D.allTimes < 100 | D.allTimes > size(self.dataSource,1)-100;
%             D.allTimes(ridx) = [];
%             D.allSessions(ridx) = [];
%             D.allAmps(ridx)  = []; 
%             D.allChans(ridx) = [];
            self.preloadedData = D;
            t2 = toc(t1); fprintf(' done. (%.3f)\n', t2);
            t1 = tic;  fprintf('Initialising waveform manager and Features...');
            self.WaveformManager = mysort.wf.MultiSessionBufferedWfManager(...
                self.dataSource, D.allTimes, D.allChans, D.allSessions, 20, 80, [], ME);
            self.Features = hdmeagui.FeatureContainer(hdmeagui.MultiSessionElectrodeFeature.empty());
            for i=1:nC
                idx = find(D.allChans==en(i));
                self.Features.featureList(i, 1) = ...
                    hdmeagui.MultiSessionElectrodeFeature(num2str(en(i)), ...
                        D.allTimes(idx), idx,...
                        D.allSessions(idx), ...
                        D.allAmps(idx), ...
                        ME, en(i));
            end
            
            t2 = toc(t1); fprintf(' done. (%.3f)\n', t2);
            self.checkSpikeSortingContainer();
        end
        %------------------------------------------------------------------
        function setDataSource(self, ds)
            assert(isa(ds, 'mysort.ds.DataSourceInterface'), 'ds must implement the DataSourceInterface!');
            self.dataSource = ds;
            self.reset();
            self.update();
        end        
        %------------------------------------------------------------------
        function CBAcceptTemplateButton(self, a, b)
            if isempty(self.SpikeSortingContainer)
                warning('You cant accept a template if no sorting is selected!');
                return
            end
            if isempty(self.SelectedWaveformManager) || isempty(self.SelectedWaveformManager.eventTimes)
                warning('No waveforms selected!');
                return
            end
%             G = mysort.mea.God(self.dataSource, self.SpikeSortingContainer, varargin)
            id = length(self.SpikeSortingContainer.unitNames) + 1;
            max_spikes_for_template_calculation = 500;
            rperm = randperm(length(self.SelectedWaveformManager.eventTimes));
            idx = rperm(1:min(length(self.SelectedWaveformManager.eventTimes), max_spikes_for_template_calculation));
            t1 = tic;  fprintf('Cutting Waveforms for Template calculation...');
            [T, ME, wfsAllSessions, nWfsPerElectrode] = self.SelectedWaveformManager.getTemplate4Idx(idx);
            t2 = toc(t1); fprintf(' done. (%.3f)\n', t2);
            
            t1 = tic;  fprintf('Initializing Template...');
            T = mysort.mea.Template(T', self.SelectedWaveformManager.wfbufferCutLeft, self.dataSource.getSamplesPerSecond(), id, ME, nWfsPerElectrode);
            t2 = toc(t1); fprintf(' done. (%.3f)\n', t2);
            es = self.SelectedWaveformManager.eventSessions;
            et = self.SelectedWaveformManager.eventTimes;
            us = unique(es);
            gdfL = cell(1, self.dataSource.getNSessions());
            for i=1:length(us)
                s = us(i);
                ets = et(es==s);
                gdfL{s} = [repmat(id, length(ets), 1) ets(:)];
            end
            self.SpikeSortingContainer.addTemplateList(gdfL, T);
            
            t1 = tic;  fprintf('Saving updated sorting...');
            self.SpikeSortingContainer.save();
            t2 = toc(t1); fprintf(' done. (%.3f)\n', t2);
            self.CBNewSorting();
            
            t1 = tic;  fprintf('Removing Spikes from detected Spikes...');
            self.removeSpikesFromFeaturesThatAreTooCloseToAlreadySortedSpikes();
            t2 = toc(t1); fprintf(' done. (%.3f)\n', t2);
            
            t1 = tic;  fprintf('Updating Featureplot...');
            self.FeaturePlot.setFeatureContainer(self.Features);  
            t2 = toc(t1); fprintf(' done. (%.3f)\n', t2);
        end
    end
end
