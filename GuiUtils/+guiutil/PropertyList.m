classdef PropertyList < handle
    properties
        % gui elements
        m_listProperties
        m_hfuncCallback
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = PropertyList(PropertyList, varargin)
            P.callback = [];
            P = mysort.util.parseInputs(P, varargin);
            self.setProperties(PropertyList, P.callback);
        end
        %------------------------------------------------------------------
        function setProperties(self, PL, cb)
            if iscell(PL);
                self.m_listProperties = guiutil.Property.empty();
                for i=1:size(PL,1)
                    self.m_listProperties(end+1) = guiutil.Property(PL{i,:}, 'changeCallback', cb);
                end
            else
                self.m_listProperties = PL;
            end
            self.m_hfuncCallback = cb;
        end
        
        %------------------------------------------------------------------
        function p = getProperty(self, pname)
            nProps = length(self.m_listProperties);
            % build strings
            for i=1:nProps
                p = self.m_listProperties(i);
                if strcmp(p.getName(), pname)
                    return
                end
            end
        end
        %------------------------------------------------------------------
        function val = getPropertyValue(self, pname)
            p = self.getProperty(pname);
            val = p.get();
        end
        %------------------------------------------------------------------
        function val = getPropertyString(self, pname)
            p = self.getProperty(pname);
            val = p.getString();
        end

      
        %------------------------------------------------------------------
        function props = getPropertyCell(self) 
            nProps = length(self.m_listProperties);
            props = {};
            for i=1:nProps
                p = self.m_listProperties(i);
                props{i,1} = p.getName();
                props{i,2} = p.get();
            end
        end
    end
end