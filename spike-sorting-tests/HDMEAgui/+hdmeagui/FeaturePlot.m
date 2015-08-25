classdef FeaturePlot < guiutil.AmplitudeSelectableAxes
    properties
       % GUI members
       featureSelectionCallback
       
       % internal members
       featureContainer
       featureNames
       featureMarkerHandle
       markedFeatureIdx 
    end
    
    methods
        %------------------------------------------------------------------
        function self = FeaturePlot(featureContainer, featureSelectionCallback, varargin)
            P.dummy = 1;
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);
            self = self@guiutil.AmplitudeSelectableAxes([], uP{:}, ...
                'clickableTicks', 1:length(featureContainer.featureList), ...
                'markerColor', 'c');
            self.featureSelectionCallback = featureSelectionCallback;
            self.setFeatureContainer(featureContainer);
            self.amplitudeSelectionCallback = @self.CBamplitudeSelection;
        end
        %------------------------------------------------------------------
        function setFeatureContainer(self, FC)
            self.featureContainer = FC;
            self.markedFeatureIdx = [];
            if isempty(FC.featureList)
                self.amplitudeSelectableTicks = [];
                self.featureNames = [];
            else
                self.amplitudeSelectableTicks = 1:length(FC.featureList);
                self.featureNames = cell(1, length(FC.featureList));
                for i=1:length(FC.featureList)
                    self.featureNames{i} = FC.featureList(i).name;
                end
            end
            self.plotFeatures();
        end
        
        %------------------------------------------------------------------
        function plotFeatures(self, FC)
            nF = self.getNFeatures();
            cla(self.Axes);
            if nF == 0
                return
            end
            set(self.Axes, 'NextPlot', 'add');
%             mysort.plot.zoomReset(self.Axes);
            t1 = tic;  fprintf('init...');
            minmax = [0 0];
            nfstotal = 0;
            if nargin == 1
                FC = self.featureContainer;
            end
            
            for i=1:nF
                f = FC.featureList(i);
                if length(f.eventFeatures)>0
                    nfstotal = nfstotal + length(f.eventFeatures);
                    [mi(i) miidx(i)] = min(f.eventFeatures);
                    [ma(i) maidx(i)] = max(f.eventFeatures);
                    minmax = [min(minmax(1), mi(i)) max(minmax(2), ma(i))];
                end
            end
            t2 = toc(t1); fprintf(' done. (%.3f)\n', t2);
            
            t1 = tic;  fprintf('building feature matrix...');
            steps = 100;
            F = zeros(steps, nF);
            range = linspace(minmax(1), minmax(2), steps);
            for i=1:nF
                f = FC.featureList(i);
                F(:, i) = histc(f.eventFeatures, range);
            end
            t2 = toc(t1); fprintf(' done. (%.3f)\n', t2);
%             if isa(f, 'hdmeagui.MultiSessionElectrodeFeature')
%                 for i=1:nF
%                     card = f.getElectrodeCardinality();
%                     if card == 1
%                         c = 'k';
%                     elseif card < f.ME.getNSessions()
%                         c = 'g';
%                     else
%                         c = 'r';
%                     end
%                     h = plot(i,0, '*', 'color', c);
%                     set(h, 'HitTest', 'off');
%                 end
%                 
%             end            
            t1 = tic;  fprintf('plotting...');
            h = imagesc(log(F+1), 'Parent', self.Axes, 'ydata', range); 
            colormap(gray)
            set(h, 'HitTest', 'off');            
            if isa(f, 'hdmeagui.MultiSessionElectrodeFeature')
                for i=1:nF
                    card = FC.featureList(i).getElectrodeCardinality();
                    if card == 1
                        c = 'w';
                        marker = '.';
                    elseif card < f.ME.getNSessions()
                        c = 'g';
                        marker = '+';
                    else
                        c = 'r';
                        marker = '*';
                    end
                    h = plot(self.Axes, i, mi(i), marker, 'color', c, 'markersize', 4, 'linewidth', 2);
                    set(h, 'HitTest', 'off');
                end
            else           
                h = plot(self.Axes, mi, '.w', 'markersize', 4, 'linewidth', 2);
                set(h, 'HitTest', 'off');            
            end
            set(self.Axes, 'xlim', [-1 nF+1]);
            set(self.Axes, 'ylim', [minmax(1) -10]);
            drawnow
            disp(minmax)
            self.plotFeatureMarker();
            t2 = toc(t1); fprintf(' done. (%.3f)\n', t2);
            
%             figure;
%             axes
%             hold on
%             FF = F(1:65,:);
%             FF(FF<2) = nan;
%             C = {};
%             for i=1:size(FF,2)
%                 plot3(1:size(FF,1),ones(1, size(FF,1))*i, log(FF(:,i)), 'color', mysort.plot.vectorColor(i));
%                 C{i} = mysort.plot.vectorColor(i);
%             end
%             mysort.plot.mc(log(FF'), 'spacer', 1, 'color', C);            
        end       
        
        %------------------------------------------------------------------
        function plotFeatureMarker(self)
            self.removeFeatureMarker();
            mfidx = self.markedFeatureIdx;
            if isempty(mfidx)
                return
            end
            self.featureMarkerHandle = plot(self.Axes, mfidx, zeros(length(mfidx),1), '*', 'color', 'm', 'markersize', 8, 'linewidth', 2);
            for i=1:length(mfidx)
                self.featureMarkerHandle(i+1) = text(mfidx(i), 0, self.featureContainer.featureList(mfidx(i)).name, 'Parent',self.Axes);
            end
            set(self.featureMarkerHandle, 'HitTest', 'off');
        end
        %------------------------------------------------------------------
        function removeFeatureMarker(self)
            if ishandle(self.featureMarkerHandle)
                try
                    delete(self.featureMarkerHandle)
                catch
                end
            end
        end
        %------------------------------------------------------------------
        function setMarkedFeatureElectrodes(self, elNumbers)
            [a fIdx c] = intersect([self.featureContainer.featureList.elnumber], elNumbers);
            self.markedFeatureIdx = fIdx;
            self.plotFeatureMarker();
        end
        
        %------------------------------------------------------------------
        function n = getNFeatures(self)
            n = length(self.featureContainer.featureList);
        end
        %------------------------------------------------------------------
        function b = CBamplitudeSelection(self, fetIdx, thr)
            if isempty(self.featureContainer.featureList) || isempty(self.featureSelectionCallback)
                b = false;
                return
            end
            b = true;
            F = self.featureContainer.featureList(fetIdx);            
            self.featureSelectionCallback(F, thr);
        end
        %------------------------------------------------------------------
        function update(self)
%             nF = self.getNFeatures();
%             for i=1:nF
%                 f = self.featureList(i);
%                 x = f.get();
% 
%             end
        end        
    end
end