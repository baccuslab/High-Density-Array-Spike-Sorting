classdef AmplitudeSelectableAxes < guiutil.SingleAxes
    properties
        amplitudeMarker
        amplitudeMarkerWidth
        amplitudeMarkerColor
        
        amplitudeSelectionCallback
        amplitudeSelectableTicks
    end
 
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = AmplitudeSelectableAxes(callback, varargin)
            P.clickableTicks = [];
            P.markerWidth = 1;
            P.markerColor = 'r';
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);
                self = self@guiutil.SingleAxes(uP{:});
            self.amplitudeSelectionCallback = callback;
            self.amplitudeMarker = [];
            self.amplitudeMarkerWidth = P.markerWidth;
            self.amplitudeMarkerColor = P.markerColor;
            self.amplitudeSelectableTicks = P.clickableTicks;
            self.makeLayout();
            set(self.Axes, 'ButtonDownFcn', @self.CBAmplitudeSelection);
            set(self.Axes, 'nextPlot', 'add');
            
        end        
        %------------------------------------------------------------------
        function CBAmplitudeSelection(self, a, b)
            point = get(gca, 'CurrentPoint');
            if isempty(self.amplitudeSelectableTicks)
                chan = round(point(1,1));                
            else
                [d midx] = min(abs(self.amplitudeSelectableTicks-point(1,1)));
                chan = self.amplitudeSelectableTicks(midx);
            end
            thr = point(1,2);
            if ~isempty(self.amplitudeSelectionCallback)
                b = self.amplitudeSelectionCallback(chan, thr);
                if ~b
                    return
                end
            end
            if ishandle(self.amplitudeMarker)
                try
                    delete(self.amplitudeMarker);
                catch
                end
            end
            set(self.Axes, 'nextplot', 'add');
            xlim = get(self.Axes, 'xlim');
            ylim = get(self.Axes, 'ylim');
            self.amplitudeMarker = ...
                plot(chan+[-.5 .5]*self.amplitudeMarkerWidth, thr+[0 0], ...
                '-', 'linewidth', 2, 'color', self.amplitudeMarkerColor);
            set(self.amplitudeMarker, 'HitTest', 'off');
            set(self.Axes, 'xlim', xlim);
            set(self.Axes, 'ylim', ylim);
        end
        %------------------------------------------------------------------
        function varargout = plot(self, varargin)
            [varargout{1:nargout}] = plot(self.Axes, varargin{:});
            set(self.Axes, 'ButtonDownFcn', @self.CBAmplitudeSelection);
        end
    end
end