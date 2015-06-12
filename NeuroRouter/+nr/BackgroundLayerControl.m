classdef BackgroundLayerControl < handle
    properties(Constant=true)
       
    end
    properties
        % gui elements
        Parent
        mainLayout
        ControlPanel        
        PropertyPanel
                
        % Views
        ChipViewHandle
        
        % callbacks
        postCallbackCheck
        
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function gui = BackgroundLayerControl(varargin)
            P.Parent = [];
            P = mysort.util.parseInputs(P, varargin, 'error');
            
            if nargin==0 || (isempty(P.Parent) &&  ~(ishandle(varargin{1}) || ...
                              isa(varargin{1}, 'uiextras.Container')))
                P.Parent = figure();
            end       
            
            gui.Parent = P.Parent;
            gui.buildLayout();
        end
        %------------------------------------------------------------------
        function buildLayout(self)
            self.mainLayout = uiextras.VBoxFlex( 'Parent', self.Parent, 'Spacing', 3 );
            self.ControlPanel = uiextras.BoxPanel( ...
                'Parent', self.mainLayout, ...
                'Title', 'Background Controls: '); 
                controlButtonGrid = uiextras.Grid( 'Parent', self.ControlPanel, 'Spacing', 1 );
            	uicontrol( 'Style', 'PushButton', ...
                	'Parent', controlButtonGrid, ...
                    'String', 'Reset', ...
                    'Callback', @(x,y) self.CBResetButton(x,y));
                set(controlButtonGrid, 'ColumnSizes', [-1 -1], 'RowSizes', repmat(20, 1, 1));
            
        end
        %------------------------------------------------------------------
        function buildPlotNeuronControl(self, parent)
        end
        
        
        %------------------------------------------------------------------
        %% CALLBACKS
        function CBResetButton(self, hObject, eventdata)
        end
    end
end

