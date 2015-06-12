classdef TMGui < guiutil.GuiElement
    properties
        % gui elements
        DataPlot
        ResidualPlot
        FilterOutputPlot
        % gui internals
        
        % internals
        dataSource
        preloadedData
        Features
    end
 
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = TMGui(varargin)
            P.dataSource = [];
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);
            self = self@guiutil.GuiElement(uP{:});
            self.dataSource = mysort.ds.checkDataSources(P.dataSource);
        end
        %------------------------------------------------------------------
        function makeLayout(self)
            p = self.Parent;            
            mainLayout = uiextras.HBox( 'Parent', p,'Padding', 0 );
                controlLayout = uiextras.VBox( 'Parent', mainLayout,'Padding', 0 ); 
                    sortingButtonBox = uiextras.BoxPanel( ...
                           'Parent', controlLayout, ...
                           'Title', 'Sorting');
                    templateBox = uiextras.BoxPanel( ...
                           'Parent', controlLayout, ...
                           'Title', 'Templates');                       
                set(controlLayout, 'Sizes', [250,-1]);   
                viewLayout = uiextras.VBoxFlex( 'Parent', mainLayout,'Padding', 0 ); 
                    self.FeaturePlot = hdmeagui.FeaturePlot([], @self.CBfeatureSelection, 'Parent', viewLayout);
                    self.SpikeWaveformPlot = mysort.plot.Spikes('Parent', viewLayout);
                    selectionDetailsBox = uiextras.HBox( 'Parent', viewLayout,...
                    'Padding', 0 );
                        Box1 = uiextras.BoxPanel( ...
                           'Parent', selectionDetailsBox, ...
                           'Title', 'Box1');
                        Box2 = uiextras.BoxPanel( ...
                           'Parent', selectionDetailsBox, ...
                           'Title', 'Box2');
                        Box3 = uiextras.BoxPanel( ...
                           'Parent', selectionDetailsBox, ...
                           'Title', 'Box3');
                    set( selectionDetailsBox, 'Sizes', [-2,-1,-1]);
                set( viewLayout, 'Sizes', [-1,-1,-2]);
             set(mainLayout, 'Sizes', [100,-1]);
        end
        %------------------------------------------------------------------
        function update(self)
            if ~self.bIsVisible
                return
            end
            disp('SortGui Update');
            self.preloadData();
%              % plot brushable plots and link data
%             h = plot(self.FeaturePlot.Axes, self.preloadedData.allChans, self.preloadedData.allAmps, ...
%                 'kx', 'HitTest', 'off', 'markersize',3);
%             axis(self.FeaturePlot.Axes, 'tight');
%             box off
% %             maxi = max(amps(idx));
% %             set(ax, 'xlim', [0 G.WF.nC+1]);
% %             set(ax, 'ylim', [0 maxi*1.01]);
% %             xlabel(ax, 'channel index');
% %             ylabel(ax, 'amplitude [std noise]');
% %             set(ax, 'color', 'none');
        end  
        %------------------------------------------------------------------
        function reset(self)
            self.preloadedData = [];         
        end  
        %------------------------------------------------------------------
        function CBfeatureSelection(self, F, thr)
            FW = hdmeagui.FeatureWaveforms(F, 20, 60);
            SP = FW.getWfs4Threshold(thr);
            disp(SP)
        end
        %------------------------------------------------------------------
        function preloadData(self)
            if ~isempty(self.preloadedData)
                return
            end
            D = self.dataSource.getSelectedSessionsDetectedSpikes();
            ME = self.dataSource.getSelectedSessionsMergedMultiElectrode();
            en = ME.getElectrodeNumbers();
            nC = ME.getNElectrodes();
            self.preloadedData = D;
            self.Features = hdmeagui.MultiSessionElectrodeFeature.empty();
            for i=1:nC
                self.Features(i, 1) = ...
                    hdmeagui.MultiSessionElectrodeFeature(num2str(en(i)), ...
                        D.allTimes(D.allChans==en(i)), ...
                        D.allSessions(D.allChans==en(i)), ...
                        D.allAmps(D.allChans==en(i)), ...
                        ME, i);
            end
            self.FeaturePlot.setFeatureList(self.Features);
            % Create Waveform manager
%             cutleft = 20;
%             tf = 60;
%             self.WF = mysort.wf.WfManagerBuffCut(self.dataSource, ...
%                 D.allTimes, D.allChans, cutleft, tf);
%             
%             % Register features
%             mima = mysort.wf.feature.MinMax(self.WF);
%             mima.needsUpdate = false;
%             mima.F = D.allAmps;
%             self.WF.registerFeature(mima);
%             
%             % Create Template Manager
%             self.T = mysort.wf.TemplateManager(H.WF);            
        end
        %------------------------------------------------------------------
        function setDataSource(self, ds)
            assert(isa(ds, 'mysort.ds.DataSourceInterface'), 'ds must implement the DataSourceInterface!');
            self.dataSource = ds;
            self.reset();
            self.update();
        end        
    end
end