classdef Property < handle
    properties
        m_strName
        m_strDisplayName
        m_type
        m_initval
        m_hfuncConstraint
        m_hfuncChangeCallback
        m_val
        m_old_val
    end
    
    methods 
        %------------------------------------------------------------------
        function self = Property(name, initval, varargin)       
            self.m_strName = name;
            P.constraint = @(x) true;
            P.initIdx = 1;
            P.changeCallback = [];
            P.displayName = [];
            P = mysort.util.parseInputs(P, varargin);
            
            self.m_strDisplayName = P.displayName;
            if isempty(P.displayName)
                self.m_strDisplayName = self.m_strName;
            end
            self.m_initval = initval;
            self.m_hfuncConstraint = P.constraint;
            if iscell(initval)
                type = 'cell';
                initval = initval{P.initIdx};
            elseif isnumeric(initval)
                type = 'numeric';
            elseif islogical(initval)
                type = 'bool';
            elseif ischar(initval)
                type = 'char';
            elseif iscell(initval)
                type = 'cell';
            else
                type = 'unknown';
            end        
            self.m_type = type;
            self.set(initval);

            self.m_hfuncChangeCallback = P.changeCallback;
        end
        
        %------------------------------------------------------------------
        function v = get(self)
            v = self.m_val;
        end
        %------------------------------------------------------------------
        function b = set(self, val)
            b = true;
            if ~isempty(self.m_hfuncConstraint)
                try
                    b = self.m_hfuncConstraint(val);
                catch ME
                    throw('CONSTRAINTFUNCTION:ERROR');
                end
                if ~b
                    return
                end
            end
            self.m_old_val = self.m_val;
            self.m_val = val;
            if ~isempty(self.m_hfuncChangeCallback)
                self.m_hfuncChangeCallback(self);
            end
        end        
        %------------------------------------------------------------------
        function v = getString(self)
            if ischar(self.m_val)
                v = self.m_val;
            elseif iscell(self.m_initval)
                v = self.m_val;
            else
                v = num2str(self.m_val);
            end
        end
        %------------------------------------------------------------------
        function str = getName(self)
            str = self.m_strName;
        end
        %------------------------------------------------------------------
        function str = getDisplayName(self)
            str = self.m_strDisplayName;
        end      
        %------------------------------------------------------------------
        function setDisplayName(self, str)
            self.m_strDisplayName = str;
        end            
    end
end