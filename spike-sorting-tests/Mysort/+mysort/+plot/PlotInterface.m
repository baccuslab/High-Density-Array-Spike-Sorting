classdef PlotInterface < mysort.util.DebuggableClass
    properties (Constant)
    end
    
    properties
    end
    
    methods (Abstract)
    end
    
    methods 
        %%% ------------------------------------------------------ 
        function self = PlotInterface(varargin)
            self = self@mysort.util.DebuggableClass(varargin{:});
            self.P.ax = [];
            self.P.figure = 1;
            self.P.fig = [];
            self.P.chunk_size = 300000; 
            self.P.datafile = [];
            suppress_warning = 1;
            self.P = mysort.util.parseInputs(self.P, 'PlotInterface', varargin, suppress_warning);
            
            self.prepareFigure();
        end
        
        %%% ------------------------------------------------------ 
        function prepareFigure(self)
            if ~isempty(self.P.fig)
                figure(self.P.fig);
                clf
            elseif self.P.figure==1 && isempty(self.P.ax)
                self.P.fig = mysort.plot.figure();
            end
        end
    end
end