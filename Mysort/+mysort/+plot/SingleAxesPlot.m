
classdef SingleAxesPlot < mysort.plot.PlotInterface
    properties (Constant)
    end
    
    properties
    end
    
    methods (Abstract)
    end
    
    methods 
        %%% ------------------------------------------------------ 
        function self = SingleAxesPlot(varargin)
            self = self@mysort.plot.PlotInterface(varargin{:});
            suppress_warning = 1;
            self.P = mysort.util.parseInputs(self.P, 'SingleAxesPlot', varargin, suppress_warning);
            
            self.prepareAxes();
        end
        
        %%% ------------------------------------------------------ 
        function prepareAxes(self)
            if ~isempty(self.P.ax)
                axes(self.P.ax(1));
                %cla reset
            else
                self.P.ax = axes();
            end
        end
    end
end