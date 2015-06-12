classdef HandleObjectTestL2  < mysort.util.HandleObjectTestL1
    % Mother class of all handle classes. Adds parameter parsing, debugging
    % output and verbose settings and some mandotory save and load
    % functions
    properties (SetAccess=private)

    end
    
    properties 
        
    end    
    
    methods (Abstract)
         
    end    
    
    methods (Access=protected)
    end
    
    %%%--------------------------------------------------------------------
    methods
        % -----------------------------------------------------------------
        function self = HandleObjectTestL2(varargin)
            [p varargin] = mysort.util.HandleObjectTestL1.checkInputs(varargin);
            addParamValue(p, 'perlich', 0, @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));
            addParamValue(p,'a',@(x) disp(x),@(x) isa(x, 'function_handle'));
            
            self@mysort.util.HandleObjectTestL1(p, 'name', 'L2', varargin{:});
        end
    end
end