classdef WaveformManagerPlot < mysort.plot.WaveformPlot
    properties
        % gui elements
        
        % gui internals
        m_nRestrict2DetectionChannels
        
        m_cMeanColor
        m_cMedianColor
        
        
        wfManager
        
        % externals
        plottedChannelNumbers
    end
    
 
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = WaveformManagerPlot(wfManager, varargin)
            P.name = 'position';
            P.restrictToDetectionChannel = [];
            P.plotMeanColor = [];
            P.plotMedianColor = [];
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);       
            uP = [uP {'waveformIDs', wfManager.eventChans}];
            self = self@mysort.plot.WaveformPlot([], [], uP{:});
            self.m_nRestrict2DetectionChannels = P.restrictToDetectionChannel;
            self.m_cMeanColor = P.plotMeanColor;
            self.m_cMedianColor = P.plotMedianColor;
            self.setWfManager(wfManager);
        end
        %------------------------------------------------------------------
        function setWfManager(self, wfManager)
            self.wfManager = wfManager;
            self.wfManager.MultiElectrode = wfManager.MultiElectrode;
            if isempty(self.wfManager)
                self.reset();
            else
                self.update();
            end
        end
        %------------------------------------------------------------------
        function reset(self)
            cla(self.Axes);
            mysort.plot.zoomReset(self.Axes);
        end
        
        %------------------------------------------------------------------
        function [idx nE] = getPlotIdx(self)
            u = self.GuiProperties.getProperty('Unit').get();
            if ~strcmp(u, 'all')
                u = str2num(u);
%                 if length(self.m_nRestrict2DetectionChannels) > 1
%                     warning('plotting more than one unit at a time not implemented! defaulting to first in selection!')
%                     self.m_nRestrict2DetectionChannels = self.m_nRestrict2DetectionChannels(1);
%                 end
%                 idx = find(self.wfManager.eventChans == self.m_nRestrict2DetectionChannels);
                idx = find(self.wfManager.eventChans == u);
            else
                idx = 1:length(self.wfManager.eventTimes);
            end          
            nE = length(idx);
            maxWf = self.GuiProperties.getProperty('plotNWfs').get();
            if nE > maxWf
                if self.GuiProperties.getProperty('plotRandWfs').get()
                    rperm = randperm(nE);
                    idx = idx(rperm(1:maxWf));
                end
                idx = idx(1:maxWf);                  
                nE = maxWf;
            end   
        end
        %------------------------------------------------------------------
        function update(self)
%             error('not implemented!');
            cla(self.Axes);
            [idx nE] = self.getPlotIdx();
            if nE == 0
                return
            end              
            h1 = self.messageGettingWaveforms();
            for wfidx=1:length(self.wfManager)
                wfsc{wfidx} = self.wfManager(wfidx).getWaveform4Idx(idx); %, 10, 30);
            end
            delete(h1);
            text(.3, .5, 'Updating Plot...', 'fontsize', 20, 'Parent', self.Axes);
            drawnow
            nC2Plot = self.GuiProperties.getProperty('plotMaxChan').get();
            for wfidx=1:length(self.wfManager)
                wfs = wfsc{wfidx};
                if iscell(wfs)
                    wfs = wfs{1};
                end

                nC = self.wfManager(wfidx).getNChannels();
                self.plottedChannelNumbers = 1:nC;
                EP = self.wfManager.MultiElectrode.electrodePositions;
                if nC2Plot < nC
                    template = mysort.wf.v2m(mean(wfs, 1), nC);
                    maxChans = max(abs(template), [], 2);
                    maxChans = sortrows([maxChans (1:nC)']);

                    plotChanIdx = zeros(1,nC);
                    plotChanIdx(maxChans(end-nC2Plot+1:end,2)) = 1;
                    self.plottedChannelNumbers = find(plotChanIdx);
                    EP = EP(plotChanIdx==1,:);
                    wfs = mysort.wf.vSubChanSel(wfs, nC, plotChanIdx==1);
                    nC = sum(plotChanIdx);
                end
                wfs = mysort.wf.v4plot(wfs, nC); 
                twfs = mysort.wf.v2t(wfs, nC);
                mysort.plot.zoomReset(self.Axes);
                self.m_nPlotMode = self.GuiProperties.getProperty('plotMode').get();
                if strcmp(self.m_nPlotMode, 'horizontal')
                    self.plotHorizontal(wfs);
                elseif strcmp(self.m_nPlotMode, 'vertical')

                elseif strcmp(self.m_nPlotMode, '2D')
                    self.plot2D(twfs, EP);
                else
                    error('unknown plot mode: %s !', self.m_nPlotMode)
                end
            end
        end
        %------------------------------------------------------------------
        function plot2D(self, wfs, EP)
            set(self.Axes, 'NextPlot', 'replace');
            mysort.plot.waveforms2D(wfs, EP, 'AxesHandle', self.Axes, 'plotMedian', 1);
            set(self.Axes, 'NextPlot', 'add');
        end        
        %------------------------------------------------------------------
        function plotHorizontal(self, wfs)        
            set(self.Axes, 'NextPlot', 'replace');
            plot(self.Axes, wfs', 'g-');
            set(self.Axes, 'NextPlot', 'add');
            if ~isempty(self.m_cMedianColor)
                plot(self.Axes, median(wfs,1), 'color', self.m_cMedianColor, 'linewidth', 3);
            end
            if ~isempty(self.m_cMeanColor)
                plot(self.Axes, mean(wfs,1)', 'color', self.m_cMeanColor, 'linewidth', 3);
            end
%             for i=1:nC
%                 text((cutLength+1)*i-cutLength/1.8, -30, num2str(ME.electrodeNumbers(i)), 'Parent', self.Axes);
%             end            
            axis(self.Axes, 'tight');
            set(self.Axes, 'xlim', [0 size(wfs, 2)+1]); 
        end
        %------------------------------------------------------------------
        function h1 = messageGettingWaveforms(self)  
            h1 = text(.3, .5, 'Getting Waveforms...', 'fontsize', 20, 'Parent', self.Axes);
                set(self.Axes, 'xlim', [0 1], 'ylim', [0 1], 'nextplot', 'replace');
                drawnow
        end
        %------------------------------------------------------------------
        function initPropertyList(self)
            initPropertyList@mysort.plot.WaveformPlot(self);
            props = {'cutLength', 35
                     'cutLeft', 10};            
            cb = @self.CBPropertyChange;
            for i=1:size(props,1)
                self.PropertyList(end+1) = guiutil.Property(props{i,:}, 'changeCallback', cb);
            end       
        end        
    end
end