
classdef spikes < mysort.plot.PlotInterface
    properties
        plotHandles
    end
    methods
        function self = spikes(spikesX, varargin)
            self = self@mysort.plot.PlotInterface(varargin{:});
            self.P.classes = ones(size(spikesX,1),1);
            self.P.IDs = [];
            self.P.linewidth = 1;
            self.P.nC = 1;
            self.P.channelIDs = [];
            self.P.stacked = false;
            self.P.linkaxes = 'xy';
            self.P.plotMean = 1;
            self.P = mysort.util.parseInputs(self.P, 'plot.spikes', varargin, 1);
            
            if isempty(spikesX)
                return
            end
%             if ndims(spikesX) == 3
%                 spikesX = mysort.wf.t2v(spikesX);
%             end
            
            if ~isempty(self.P.IDs)
                self.P.classes = self.P.IDs;
            end
            if self.P.nC == 1 && ~isempty(self.P.channelIDs)
                self.P.nC = length(self.P.channelIDs);
            end
            % Shortcuts
            cls = unique(self.P.classes);
            nClasses = length(cls);
            nPlots = nClasses;
            if ~self.P.stacked && length(self.P.ax) ~= nClasses
                self.P.ax = mysort.plot.subplots([nClasses 1],'spacerY',.02);
            elseif self.P.stacked
                if length(self.P.ax) ~= 1
                    self.P.ax = mysort.plot.subplots(1);
                end
                nPlots = 1;
            end
            
            if max(size(self.P.classes))==1
                self.P.classes = repmat(self.P.classes, size(spikesX,1),1);
            end
            
           
            self.plotHandles = mysort.plot.ConcatenatedWaveform.empty(nPlots,0);
            noxtick = 1;
            for i=1:nClasses
                myId = cls(i);
                idx = self.P.classes==myId;
                if i==nClasses
                    noxtick = 0;
                end
                phidx = i;
                if nPlots == 1
                    phidx = 1;
                end
                if i == 1
                    channelIDs = self.P.channelIDs;
                else
                    channelIDs = [];
                end
                self.plotHandles(phidx) = mysort.plot.ConcatenatedWaveform(...
                    spikesX(idx,:),'Id', myId, 'ax', self.P.ax(phidx),...
                    'NoXTickLabel',noxtick, 'nC', self.P.nC, 'linewidth', self.P.linewidth,...
                    'colorWithId', true, 'channelIDs', channelIDs, 'plotMean', self.P.plotMean);
                hold on
            end
            if ~isempty(self.P.linkaxes) & self.P.linkaxes
%                 linkaxes(self.P.plotHandles, self.P.linkaxes);
            end
        end
    end
end