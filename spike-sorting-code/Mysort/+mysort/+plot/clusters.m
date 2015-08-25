
classdef clusters < mysort.plot.PlotInterface
    properties
        plotHandles
    end
    methods
        function self = clusters(spikesX, classes, C, varargin)
            self = self@mysort.plot.PlotInterface(varargin{:});
            warning('Dont use this file, use clusterPCA instead!');
            self.P.markersize = 12;
            self.P.templates = [];
            self.P.IDs = [];
            self.P.iU = [];
            self.P = mysort.util.parseInputs(self.P, 'plot.clusters', varargin, 1);
            
            if nargin < 2 || isempty(classes)
                classes = ones(size(spikesX,1),1);
            end
            if nargin < 3
                C = [];
            end
            
            % Shortcuts
            if isempty(self.P.IDs)
                self.P.IDs = unique(classes);
            end
            nClasses = length(self.P.IDs);
            if isempty(self.P.iU) && ~isempty(C)
                U = chol(C);
                self.P.iU = inv(U);
                pwSpikesX = spikesX*self.P.iU;
            else
                pwSpikesX = spikesX;
            end
            
            
            if length(self.P.ax) ~= nClasses
                self.P.ax = mysort.plot.subplots(nClasses,'spacerY',.02);
                mysort.plot.figureName('IntraClusterPCA');
            end
            for i=1:nClasses
                myId = self.P.IDs(i);
                idx = classes==myId;
                [fetX pcs] = mysort.util.dimReductionPCA(pwSpikesX(idx,:), 2);
                set(self.P.ax(i), 'nextPlot', 'replace')
                plot(self.P.ax(i), fetX(:,1), fetX(:,2), '.', 'color', mysort.plot.vectorColor(myId), 'markersize', self.P.markersize);
                set(self.P.ax(i), 'nextPlot', 'add')
%                 if ~isempty(self.P.templates)
%                     myTemp = self.P.templates(i,:)*iU;
%                     myTF = myTemp*pcs(1:2,:)';
%                     mysort.plot.circle(myTF, 3, 'color', [0 0 0], 'linewidth', 2);
%                     plot(myTF(1), myTF(2), 'xk', 'markerSize', 10, 'linewidth', 2); 
%                 else
%                     mysort.plot.circle([0 0], 3, 'color', [0 0 0], 'linewidth', 2);
%                 end
                mysort.plot.circle(self.P.ax(i), [0 0], 3, 'color', [0 0 0], 'linewidth', 2);
                plot(self.P.ax(i), 0, 0, 'xk', 'markerSize', 10, 'linewidth', 1.5); 
            end
        end
    end
end 