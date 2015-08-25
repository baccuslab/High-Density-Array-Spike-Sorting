classdef MultiSessionWaveformManagerPlot < mysort.plot.WaveformManagerPlot
    properties
        % gui elements
        
        % gui internals
        T
        ME
        waveforms
    end
    
 
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = MultiSessionWaveformManagerPlot(wfManager, varargin)
            P.dummy = [];
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);       
            self = self@mysort.plot.WaveformManagerPlot(wfManager, uP{:});
        end
       
        %------------------------------------------------------------------
        function update(self)
            cla(self.Axes);
            if ~self.bIsVisible || isempty(self.wfManager)
                return
            end               
            [idx nE] = self.getPlotIdx();
            if nE == 0
                return
            end
            self.plotWaveforms4Idx(idx);
        end
        %------------------------------------------------------------------
        function plotWaveforms4Idx(self, idx, cutLeft, cutLength, electrodeNumbers)
            h1 = self.messageGettingWaveforms();
            if nargin <3 || isempty(cutLeft)
                cutLeft = self.wfManager.wfbufferCutLeft;
            end
            if nargin <4 || isempty(cutLength)
                cutLength = self.wfManager.wfbufferCutLength;
            end
            if nargin <5 || isempty(electrodeNumbers)
                gEN = self.wfManager.MultiElectrode.electrodeNumbers;
            else
                gEN = electrodeNumbers;
            end
            [T, ME, wfsAllSessions] = self.wfManager.getTemplate4Idx(idx, cutLeft, cutLength, gEN);
            self.T = T;
            self.ME = ME;
            self.waveforms = wfsAllSessions;
            nC = ME.getNElectrodes();
            self.plottedChannelNumbers = ME.electrodeNumbers;
            plotMaxChan = self.GuiProperties.getProperty('plotMaxChan').get();
            if plotMaxChan < nC
                maxChans = max(abs(T), [], 2);
                maxChans(isnan(maxChans)) = -1;
                maxChans = sortrows([maxChans (1:nC)']);
                
                plotChanIdx = sort(maxChans(end-plotMaxChan+1:end,2));
  
                elNumbers = ME.electrodeNumbers(plotChanIdx);
                self.plottedChannelNumbers = elNumbers;
                nC = length(elNumbers);
                [T, ME, wfsAllSessions] = self.wfManager.getTemplate4Idx(idx, cutLeft, cutLength, elNumbers);
            end
            
            delete(h1);
            text(.3, .5, 'Updating Plot...', 'fontsize', 20, 'Parent', self.Axes);
            drawnow
            mysort.plot.zoomReset(self.Axes);
            set(self.Axes, 'NextPlot', 'replace');
            
            plotMode = self.GuiProperties.getProperty('plotMode').get();
            if strcmp(plotMode, 'horizontal')
                plotHorizontal();
            elseif strcmp(plotMode, 'vertical')
                plotVertical();
            elseif strcmp(plotMode, '2D')
                plot2d();
            else
                error('unknown plot mode: %s !', self.m_nPlotMode)
            end
            %------------------------------------------------------------------
            function plotHorizontal()
                for i=1:length(wfsAllSessions)
                    wfs = wfsAllSessions{i};
                    if isempty(wfs)
                        continue
                    end
                    myNC = ME.MultiElectrodeList(i).getNElectrodes();
                    wfs = mysort.wf.v4plot(wfs, myNC);
                    rangeidx = mysort.wf.vSubChannelIdx(cutLength+1, nC, ME.SingleSessionElectrodeIndices{i});

                    if i==2
                        LineStyle = ':';
                    else
                        LineStyle = '-';
                    end

                    plot(self.Axes, rangeidx, wfs', 'color', self.waveformColor, 'LineStyle', LineStyle);
                    set(self.Axes, 'NextPlot', 'add');
                    if ~isempty(self.m_cMedianColor)
                        plot(self.Axes, rangeidx, median(wfs,1), 'color', self.m_cMedianColor, 'linewidth', 3, 'LineStyle', LineStyle);
                    end
                    if ~isempty(self.m_cMeanColor)
                        plot(self.Axes, rangeidx, mean(wfs,1)', 'color', self.m_cMeanColor, 'linewidth', 3, 'LineStyle', LineStyle);
                    end
                end
                if self.GuiProperties.getProperty('plotTemplate').get()
                    plot(self.Axes, mysort.wf.m2v([T nan(nC,1)]), 'color', self.templateColor, 'linewidth', 3);
                end
                
                if self.GuiProperties.getProperty('plotElNumbers').get()
                    for i=1:nC
                        text((cutLength+1)*i-cutLength/1.8, -30, num2str(ME.electrodeNumbers(i)), 'Parent', self.Axes);
                    end
                end
                axis(self.Axes, 'tight');
    %             set(self.Axes, 'xlim', [0 size(wfs, 2)+1]);
            end
            %------------------------------------------------------------------
            function plotVertical()
               error('not implemented yet');
            end
            %------------------------------------------------------------------
            function plot2d()
                for i=1:length(wfsAllSessions)
                    wfs = wfsAllSessions{i};
                    if isempty(wfs)
                        continue
                    end
                    nC = ME.MultiElectrodeList(i).getNElectrodes();
                    EP = ME.MultiElectrodeList(i).electrodePositions;
                    self.plotWaveforms(mysort.wf.v2t(wfs, nC), ...
                        EP, 'plotTemplate', false);                                    
%                     myNC = ME.MultiElectrodeList(i).getNElectrodes();
%                     wfs = mysort.wf.v4plot(wfs, myNC);
%                     rangeidx = mysort.wf.vSubChannelIdx(cutLength+1, nC, ME.SingleSessionElectrodeIndices{i});
% 
%                     mysort.
                end
                self.plotWaveforms([], ME.electrodePositions, 'template', T');
                Tf = size(T,1);
                if strcmp(self.GuiProperties.getProperty('markElCardinality').get(), '*')
                    idx = ME.electrodeCardinality==1;
                    plot(self.Axes, ME.electrodePositions(idx,1), ME.electrodePositions(idx,2),...
                        '*', 'color', [.8 .8 .8], 'markersize', 14, 'linewidth', 2);
                    if length(ME.MultiElectrodeList)>1
                        idx = ME.electrodeCardinality > 1 & ME.electrodeCardinality < length(ME.MultiElectrodeList);
                        plot(self.Axes, ME.electrodePositions(idx,1), ME.electrodePositions(idx,2), '*', 'color', 'g', 'markersize', 14, 'linewidth', 2);                    
                    
                        idx = ME.electrodeCardinality == length(ME.MultiElectrodeList);
                        plot(self.Axes, ME.electrodePositions(idx,1), ME.electrodePositions(idx,2), '*', 'color', 'r', 'markersize', 14, 'linewidth', 2);                    
                    end
                end
%                 plot(self.Axes, mysort.wf.m2v([T nan(nC,1)]), 'color', 'b', 'linewidth', 3);
%                 for i=1:nC
%                     text((cutLength+1)*i-cutLength/1.8, -30, num2str(ME.electrodeNumbers(i)), 'Parent', self.Axes);
%                 end
            end  
        end
        %------------------------------------------------------------------
        function initPropertyList(self)
            initPropertyList@mysort.plot.WaveformManagerPlot(self);
            props = {'markElCardinality', {'off', '*'}};            
            cb = @self.CBPropertyChange;
            for i=1:size(props,1)
                self.PropertyList(end+1) = guiutil.Property(props{i,:}, 'changeCallback', cb);
            end                 
        end
    end
end