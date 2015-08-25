
classdef DebuggableClass < handle
    properties (Constant)
        LEVEL_EXTREME = 15;
        LEVEL_PROCESSING = 10;
        LEVEL_CHUNKING = 8;
        LEVEL_FLOW = 4;
        LEVEL_WARNING = 2;    
    end
    properties
        was_restored
        P
    end
    
    methods (Abstract)

    end
    
    methods
        %%% ------------------------------------------------------
        function self = DebuggableClass(varargin)
            self.was_restored = false;
            if nargin == 2 && strcmp(varargin{1}, 'RESTORE_FROM_STRUCT')
                self.restoreFromStruct(varargin{2});
                self.was_restored = true;
                return
            end
            self.P.verbose = self.LEVEL_PROCESSING;
            self.P.log_function = @(x) disp(x);
            self.P.debug_name = '';
            suppress_warning = 1;
            self.P = mysort.util.parseInputs(self.P, 'DebuggableClass', varargin, suppress_warning);         
        end
             
        %%% ------------------------------------------------------ 
        function set(self, varargin)
            self.P = mysort.util.parseInputs(self.P, 'DebuggableClass.set', varargin);
        end
        
        %%% ------------------------------------------------------ 
        function restoreFromStruct(self, S)
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
        function S = getStruct(self)
            fields = fieldnames(self);
            for f=1:length(fields)
                S.(fields{f}) = self.(fields{f});
            end
        end        
        %%% ------------------------------------------------------
        %%% -----------------HELPER FUNCTIONS---------------------
        %%% ------------------------------------------------------
        function debugout(self, str, level)
            if ~isempty(self.P.debug_name)
                str = [self.P.debug_name ': ' str];
            end
            if self.P.verbose >= level
                self.P.log_function(str);
            end
        end 
    end
end