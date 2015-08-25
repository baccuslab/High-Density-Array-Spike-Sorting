classdef TemplateMatchingGui < guiutil.GuiElement
    properties
        % gui elements
        Axes
        RunTMButton
        NoButton
        % gui internals
        
        % internals
        dataSource
        multiSessionDataSource
        BOTM
        SourceSpikeSortingContainer
        TargetSpikeSortingContainer
        
        % callbacks 
        CBNewSorting        
    end
 
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = TemplateMatchingGui(varargin)
            P.dataSource = [];
            P.newSortingCallback = [];
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);
            self = self@guiutil.GuiElement(uP{:});  
            self.makeLayout();
            self.setDataSource(P.dataSource);
            self.CBNewSorting = P.newSortingCallback;
            
        end
        %------------------------------------------------------------------
        function makeLayout(self)
            p = self.Parent;                 
            mainLayout = uiextras.VBoxFlex( 'Parent', p, 'Spacing', 3 );
                controlButtonBox = uiextras.BoxPanel( ...
                               'Parent', mainLayout, ...
                               'Title', 'Controls');
                   controlButtonLayout = uiextras.HBox( 'Parent', controlButtonBox, 'Padding', 0 ); 
                       self.RunTMButton = uicontrol( 'Style', 'PushButton', ...
                               'Parent', controlButtonLayout, ...
                               'String', 'Run BOTM', ...
                               'Callback', @self.CBRunTMButton); 
                       self.NoButton = uicontrol( 'Style', 'PushButton', ...
                               'Parent', controlButtonLayout, ...
                               'String', 'Dummy', ...
                               'Callback', @self.CBRunTMButton); 
                       uiextras.Empty( 'Parent', controlButtonLayout )
                   set(controlButtonLayout, 'Sizes', [100, 100, -1]);
                TMbox = uiextras.BoxPanel( ...
                               'Parent', mainLayout, ...
                               'Title', 'Template Matching');
                    self.Axes = mysort.plot.SliderDataAxes({}, 'Parent', TMbox,...
                        'plotSortings', 1, 'channelSpacers', [30 0], 'zeroLineStyle', {[], ':k'});
            set( mainLayout, 'Sizes', [50, -1]);                    
        end
        %------------------------------------------------------------------
        function setDataSource(self, ds)
            if isempty(ds)
                self.dataSource = [];
                self.reset();
                return
            end
%             if isa(ds, 'mysort.ds.MultiSessionInterface')
%                 self.multiSessionDataSource = ds;
%                 ds = ds.getActiveSession();
%             end
            assert(isa(ds, 'mysort.ds.PreProcessedDataSourceInterface'), 'ds must implement the PreProcessedDataSourceInterface!');
            self.dataSource = ds;
            self.update();
        end
        %------------------------------------------------------------------
        function setSelectedSortingIndex(self, v)
%             self.SourceSpikeSortingContainer = self.multiSessionDataSource.SpikeSortingContainers{v};
%             sidx = self.multiSessionDataSource.getActiveSessionIdx();
            if ~self.bIsVisible
                return
            end
            if isempty(self.dataSource)
                warning('no datasource set in TemplateMatchinGui!');
                return
            end
            if ~isempty(self.Axes.Axes)
                for i=1:length(self.Axes.Axes)
                    cla(self.Axes.Axes(i));
                end
                set(self.Axes.Axes(1), 'xlim', [0 1], 'ylim', [0 1]);
                text(.3, .5, 'Loading...', 'Parent',self.Axes.Axes(1), 'fontsize', 20);
                drawnow
            end
            ds = self.dataSource.getActiveSession();
            self.dataSource.setActiveSpikeSortingIdx(v);
            S = self.dataSource.getActiveSpikeSorting();
            
            T = S.getTemplateWaveforms();
%             T = mysort.wf.pruneTemplates(T, ...
%                 'maxChannels', 50,...
%                 'minChannels', 2, ...
%                 'absThreshold', 2,...
%                 'setInvalidChannelsTo', 0);            
            T = mysort.wf.pruneTemplates(T, ...
                'maxChannels', 30,...
                'minChannels', 0, ...
                'absThreshold', 15,...
                'setInvalidChannelsTo', 0); 
%             self.dataSource.SpikeSortingContainers
%             S = self.SourceSpikeSortingContainer.getSingleSessionSorting();
            disp('Init BOTM...');
            Tf = size(T,1);
            Tfmax = 30;
            if Tf>Tfmax
                T = T(1:Tfmax,:,:);
                Tf = Tfmax;
            end
            Covest = ds.getCovest();
%             figure; imagesc(Covest.CCol)
%             figure; mysort.plot.waveforms2D(T, ds.MultiElectrode.electrodePositions, 'plotArgs', {'-', 'linewidth', 2});
            
            botm = mysort.sorters.BOTM(Covest, Tf, mysort.wf.t2v(T), ...
                'spikePrior', .01,...
                'max_num_of_channels_per_template', 30, ...
                'adaptOnInit', true);
            disp('BOTM done. Setting Plot...'); %'upsample', 1
            botm.DH = ds;
%             botm.adapt();
%             self.dataSource.SpikeSortingContainers = {self.SourceSpikeSortingContainer};
            self.Axes.setDataSource({self.dataSource, botm});
%             T = self.SourceSpikeSortingContainer;
            self.SourceSpikeSortingContainer = self.dataSource.SpikeSortingContainers{v};
            self.BOTM = botm;
            disp('All done.');
        end
        %------------------------------------------------------------------
        function CBRunTMButton(self, a, b)
            S = self.SourceSpikeSortingContainer;
            if isempty(S)
                warning('No source spike sorting selected!')
                return
            end
            prompt = {'Enter name for new spike sorting:'};
            dlg_title = 'Input for spike sorting name';
            num_lines = 1;
            def = {'no_name'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            if isempty(answer)     
                return
            end
            sessionIdx = self.dataSource.getSelectedSessionIdx();
            self.TargetSpikeSortingContainer = self.dataSource.createNewSpikeSortingContainer(answer{1});
            TargetSSC = mysort.sorters.runBOTMOnSpikeSortingContainer(self.dataSource,...
                sessionIdx, self.TargetSpikeSortingContainer, 30);
            if ~isempty(self.CBNewSorting)
                self.CBNewSorting(); %TargetSSC
            end            
        end
        %------------------------------------------------------------------
        function update(self)
            if ~self.bIsVisible
                return
            end
            disp('TemplateMatchingGui Update');
            if isempty(self.dataSource)
                self.reset();
                return
            end
        end
        %------------------------------------------------------------------
        function reset(self)
            self.dataSource = [];
            self.Axes.reset();
        end             
    end
end