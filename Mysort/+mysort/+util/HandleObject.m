classdef HandleObject  < handle
    % Mother class of all handle classes. Adds parameter parsing, debugging
    % output and verbose settings and some mandotory save and load
    % functions
    properties (SetAccess=private)
        inputParserObject
        p
    end
    
    properties 
        
    end    
    
    methods (Abstract)
         
    end    
    
    methods (Access=protected)
    end
    
    methods (Static, Access=protected)
        % -----------------------------------------------------------------
        function [p v] = checkInputs(v)
            if ~isempty(v) > 0 && isa(v{1}, 'inputParser')
                p = v{1};
                v(1) = [];
            else
                p = inputParser;
            end
        end
    end
    
    
    %%%--------------------------------------------------------------------
    methods
        % -----------------------------------------------------------------
        function self = HandleObject(varargin)
            [ip varargin] = mysort.util.HandleObject.checkInputs(varargin);
            addParamValue(ip, 'name', '', @(x) ischar(x));            
            addParamValue(ip, 'verbose', 0, @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));
            addParamValue(ip,'stdout_function',@(x) disp(x),@(x) isa(x, 'function_handle'));
            ip.parse(varargin{:});
            self.inputParserObject = ip;
            self.p = ip.Results;
        end 
        
        %%% ------------------------------------------------------ 
        function fromStruct(self, S)
            fields = fieldnames(S);
            for f=1:length(fields)
                try
                    self.(fields{f}) = S.(fields{f});
                catch ME
                    %warning(sprintf('Could not set property, probably constant: %s\n', fields{f}));
                end
            end
        end
        
        %%% ------------------------------------------------------ 
        function S = toStruct(self)
            fields = fieldnames(self);
            for f = 1:length(fields)
                S.(fields{f}) = self.(fields{f});
            end
        end        

        %%% ------------------------------------------------------
        function stdout(self, str, level)
            if nargin < 3
                level = 1;
            end
            if self.p.verbose < level
                return
            end            
            if ~isempty(self.p.name)
                str = [self.p.name ': ' str];
            end
            self.p.stdout_function(str);
        end    
    end
end