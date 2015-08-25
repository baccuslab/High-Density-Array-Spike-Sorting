classdef SpikeGui < guiutil.GuiElement
    properties
        % gui elements
        Axes
        
        % gui internals
        
        % internals
        dataSource
    end
 
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = SpikeGui(varargin)
            P.dataSource = [];
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);
            self = self@guiutil.GuiElement(uP{:});
            self.makeLayout();
            self.setDataSource(P.dataSource);
        end
        %------------------------------------------------------------------
        function makeLayout(self)
            p = self.Parent;           
            mainLayout = uiextras.HBoxFlex( 'Parent', p, 'Spacing', 3 );
            self.Axes = axes('Parent', mainLayout);
        end
        %------------------------------------------------------------------
        function setDataSource(self, ds)
            if isempty(ds)
                self.dataSource = [];
                self.reset();
                return
            end
%             if isa(ds, 'mysort.ds.MultiSessionInterface')
%                 ds = ds.getActiveSession();
%             end
            assert(isa(ds, 'mysort.ds.PreProcessedDataSourceInterface'), 'ds must implement the PreProcessedDataSourceInterface!');
            self.dataSource = ds;
            self.update();
        end
        %------------------------------------------------------------------
        function update(self)
            if ~self.bIsVisible
                return
            end
            if isempty(self.dataSource)
                self.reset();
                return
            end
            if ~self.dataSource.isPreprocessed()
%                 warning('The datasource was not preprocessed yet. Do that before using this plot.');
                self.reset();
                text(.2, .5, 'Not preprocessed yet', 'Parent', self.Axes, 'fontsize', 20);
                return
            end
            
            ME = self.dataSource.getSelectedSessionsMergedMultiElectrode();
            nSess = self.dataSource.getNSelectedSessions();
            nEl = ME.getNElectrodes();
            allsession_el_nr = ME.getElectrodeNumbers(); 
            allsession_el_card = ME.getElectrodeCardinality();
            oneEl   = allsession_el_card == 1;
            moreEl  = allsession_el_card > 1 & allsession_el_card < nSess;
            allEl   = allsession_el_card == nSess & nSess > 1;
            
            D = self.dataSource.getSelectedSessionsDetectedSpikes();        
            offsets = self.dataSource.getSelectedSessionsTimeOffsets();
            timeMultiplier = 1000/self.dataSource.getSamplesPerSecond();
            cla(self.Axes);
            D.allTimes = D.allTimes + offsets(D.allSessions);
            plot(self.Axes, D.allTimes*timeMultiplier, D.allChans, '.', 'color', [.5 .5 .5]);  
            set(self.Axes, 'NextPlot', 'add');
            for i=1:nEl
                if oneEl(i)
%                     c = [.5 .5 .5];
                    continue
                elseif moreEl(i)
                    c = [0 1 0];
                else
                    c = [1 0 0];
                end
                plot(self.Axes, -1, allsession_el_nr(i), '.', 'color', c, 'markersize', 14)
            end               
            
            xlabel(self.Axes, 'time in ms');
            ylabel(self.Axes, 'channel index');
            
            axis(self.Axes, 'tight')
        end
        %------------------------------------------------------------------
        function reset(self)
            self.dataSource = [];
            cla(self.Axes);
        end             
    end
end