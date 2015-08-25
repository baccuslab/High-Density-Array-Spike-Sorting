classdef HandleObjectTestL1  < mysort.util.HandleObject
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
        function self = HandleObjectTestL1(varargin)
            [p varargin] = mysort.util.HandleObject.checkInputs(varargin);
            addParamValue(p, 'bla', 0, @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));
            addParamValue(p,'peter',@(x) disp(x),@(x) isa(x, 'function_handle'));
            self@mysort.util.HandleObject(p, 'name', 'L1', varargin{:});
        end
    end
end