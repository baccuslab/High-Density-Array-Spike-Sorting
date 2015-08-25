classdef SingleAxes < guiutil.GuiElement
    properties
        % gui elements
        Axes
    end
    
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = SingleAxes(varargin)
            P.dummy = 1;
            P.AxesHandle = [];
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);
            self = self@guiutil.GuiElement(uP{:});
            self.Axes = P.AxesHandle;
        end
        %------------------------------------------------------------------
        function makeLayout(self)
            if isempty(self.Axes)
                p = self.Parent;
                self.Axes = axes('Parent', p);
            end
        end
        %------------------------------------------------------------------
        function update(self)
        end
        %------------------------------------------------------------------
        function varargout = plot(self, varargin)
            [varargout{1:nargout}] = plot(self.Axes, varargin{:});
        end
        %------------------------------------------------------------------
        function varargout = get(self, varargin)
            [varargout{1:nargout}] = get(self.Axes, varargin{:});
        end
        %------------------------------------------------------------------
        function varargout = set(self, varargin)
            [varargout{1:nargout}] = set(self.Axes, varargin{:});
        end
    end
end