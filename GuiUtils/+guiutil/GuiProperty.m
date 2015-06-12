classdef GuiProperty < handle
    properties
        m_hguiParent
        m_hguiDisplayString
        m_Property
        m_type
        m_bBoolWithOwnString
        m_hGuiHandle
        m_hfuncVal2GuiVal
        m_hfuncGuiVal2Val
    end
    
    methods 
        %------------------------------------------------------------------
        function self = GuiProperty(Parent, Property, varargin)  
            self.m_hguiParent = Parent;
            self.m_Property = Property;
            initval = Property.m_initval;
            P.displayString = [];
            P.boolWithOwnString = 'true';
            if isnumeric(initval)
                type = 'numeric';
                P.constraint = @(x) isnumeric(x) && ~isempty(x);
                P.val2guival = @num2str;
                P.guival2val = @str2num;
            elseif islogical(initval)
                type = 'bool';
                P.constraint = @(x) islogical(x) && ~isempty(x);
                P.val2guival = @(x) x;
                P.guival2val = @(x) x==1;
            elseif ischar(initval)
                type = 'char';
                P.constraint = @(x) ischar(x) && ~isempty(x);
                P.val2guival = @(x) x;
                P.guival2val = @(x) x;
            elseif iscell(initval)
                type = 'cell';
                P.constraint = @(x) ~isempty(mysort.util.findCellIdx(initval, x));
                P.val2guival = @(x) mysort.util.findCellIdx(initval, x);
                P.guival2val = @(x) initval{x};
            else
                error('type must be one of those: numeric, char, bool, cell');
            end
            P = mysort.util.parseInputs(P, varargin);
            self.m_hguiDisplayString = P.displayString;
            self.m_Property.m_hfuncConstraint = P.constraint;
            self.m_type = type;
            self.m_bBoolWithOwnString = P.boolWithOwnString;
            self.m_hfuncVal2GuiVal = P.val2guival;
            self.m_hfuncGuiVal2Val = P.guival2val;
            self.createGuiHandle();
        end
        %------------------------------------------------------------------
        function createGuiHandle(self)
            p = self.m_hguiParent;
            if strcmp(self.m_type, 'bool')
                if self.m_bBoolWithOwnString
                    self.m_hGuiHandle = uicontrol('Parent', p,...
                        'Style', 'checkbox', ...
                        'String', self.m_Property.m_strDisplayName,...
                        'value', self.get(),...
                        'callback', @self.callback);
                else
                    self.m_hGuiHandle = uicontrol('Parent', p,...
                        'Style', 'checkbox', ...
                        'value', self.get(),...
                        'callback', @self.callback);                   
                end
            elseif strcmp(self.m_type, 'cell')
                self.m_hGuiHandle = uicontrol('Parent', p,...
                    'Style', 'popupmenu', ...
                    'String', self.m_Property.m_initval(),...
                    'callback', @self.callback);
            else
                self.m_hGuiHandle = uicontrol('Parent', p,...
                    'Style', 'edit', ...
                    'String', self.getString(),...
                    'callback', @self.callback);
            end
        end
        %------------------------------------------------------------------
        function callback(self, hObject, eventdata)
            if strcmp(self.m_type, 'numeric') || strcmp(self.m_type, 'char')
                guiprop = 'String';
            elseif strcmp(self.m_type, 'cell')
                guiprop = 'Value';              
            else
                guiprop = 'Value';
            end
            guival = get(hObject, guiprop);
            val = self.m_hfuncGuiVal2Val(guival);
            
            if ~self.m_Property.set(val)
                warning('conversion of property is not ok, value changed back.') 
                set(hObject, guiprop, self.m_hfuncVal2GuiVal(self.m_Property.get()));
            end
        end
        
        %------------------------------------------------------------------
        function v = get(self)
            v = self.m_Property.get();
        end
        %------------------------------------------------------------------
        function set(self, v)
            if ~self.m_Property.set(v)
                warning('conversion of property is not ok, value was not set.') 
                return
            end
            if strcmp(self.m_type, 'numeric') || strcmp(self.m_type, 'char')
                guiprop = 'String';
            elseif strcmp(self.m_type, 'cell')
                guiprop = 'Value';
                v = mysort.util.findCellIdx(get(self.m_hGuiHandle, 'String'), v);
            else
                guiprop = 'Value';
            end
            set(self.m_hGuiHandle, guiprop, v);            
        end
        
        %------------------------------------------------------------------
        function str = getString(self)
            str = self.m_Property.getString();
        end
        %------------------------------------------------------------------
        function str = getName(self)
            str = self.m_Property.getName();
        end
        %------------------------------------------------------------------
        function str = getDisplayName(self)
            str = self.m_Property.getDisplayName();
        end    
        %------------------------------------------------------------------
        function str = setDisplayName(self, str)
            self.m_Property.setDisplayName(str);
            if ishandle(self.m_hguiDisplayString)
                set(self.m_hguiDisplayString, 'String', str);
            end
        end  
    end
end