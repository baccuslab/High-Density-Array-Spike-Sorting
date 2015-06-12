classdef GuiProperties < guiutil.PropertyList
    properties
        % gui elements
        m_hguiParent
        m_hguiPropertyGrid
        LoadButton
        SaveButton
        m_strLabelArrangement
        m_strArrangement
        % vars
        hasLoadSaveButtons
        m_listGuiStrings
        m_listGuiProperties
        m_nStringLayoutWidth
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = GuiProperties(parentHandle, PropertyList, varargin)
            P.loadSaveButton = 1;
            P.arrangement = 'vertical';
            P.labelArrangement = 'horizontal';
            P.stringLayoutWidth = -2;
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);
            self = self@guiutil.PropertyList(PropertyList, uP{:});
            
            self.m_hguiParent = parentHandle;            
            self.m_listGuiProperties = [];
            self.m_nStringLayoutWidth = P.stringLayoutWidth;
            self.m_strLabelArrangement = P.labelArrangement;
            self.m_strArrangement = P.arrangement;
            self.hasLoadSaveButtons = P.loadSaveButton;
            self.buildGui();
        end
 
        %------------------------------------------------------------------
        function buildGui(self)
            if isempty(self.m_hguiPropertyGrid)
                if self.hasLoadSaveButtons
                    outerLayout =  uiextras.VBox('Parent', self.m_hguiParent, 'Spacing', 3 );
                    
                    layout = uiextras.HBox('Parent', outerLayout, 'Spacing', 3 );
                    self.LoadButton = uicontrol( 'Style', 'PushButton', ...
                       'Parent', layout, ...
                       'String', 'Load', ...
                       'Callback', @self.CBLoadButton);
                    self.SaveButton = uicontrol( 'Style', 'PushButton', ...
                       'Parent', layout, ...
                       'String', 'Save', ...
                       'Callback', @self.CBSaveButton);
                   self.m_hguiPropertyGrid = uiextras.Grid( 'Parent', outerLayout, 'Spacing', 1 );
                   set(outerLayout, 'Sizes', [20, -1]);
                else
                    self.m_hguiPropertyGrid = uiextras.Grid( 'Parent', self.m_hguiParent, 'Spacing', 1 );
                end
            end
            self.buildGuiProps();
        end
        
        %------------------------------------------------------------------
        function p = buildGuiProps(self)
            PL = self.m_listProperties;
            nProps = length(PL);
            self.m_listGuiProperties = guiutil.GuiProperty.empty();
            layout = self.m_hguiPropertyGrid;

            if strcmp(self.m_strLabelArrangement, 'horizontal')
                % build strings
                for i=1:nProps
                    p = PL(i);
                    self.m_listGuiStrings(i) = uicontrol('Parent', layout, 'Style', 'text',...
                        'String', p.getDisplayName());
                end

                % build edit fields
                for i=1:nProps
                    p = PL(i);
                    self.m_listGuiProperties(end+1) = guiutil.GuiProperty(layout, p, ...
                        'boolWithOwnString', false, 'displayString', self.m_listGuiStrings(i));
                end
                set(layout, 'ColumnSizes', [self.m_nStringLayoutWidth -1], 'RowSizes', repmat(20, nProps, 1));
            else
                % build strings and edit fields one after the other
                rowsizes = [];
                for i=1:nProps
                    p = PL(i);
                    if ~strcmp(p.m_type, 'bool')
                        self.m_listGuiStrings(i) = uicontrol('Parent', layout, 'Style', 'text',...
                            'String', p.getDisplayName());
                        rowsizes(end+1) = 15;
                    end
                    self.m_listGuiProperties(end+1) = guiutil.GuiProperty(layout, p,...
                        'boolWithOwnString', true);
                    rowsizes(end+1) = 22;
                end
                set(layout, 'ColumnSizes', [-1], 'RowSizes', rowsizes);                
            end
        end
%         %------------------------------------------------------------------
%         function self = setProperties(self, PL, cb)
%             setProperties@guiutil.PropertyList(self, PL, cb);
%             for i=1:length(self.m_listProperties)
%                 self.m_listGuiProperties(i).m_Property = self.m_listProperties(i);
%             end
%         end
        %------------------------------------------------------------------
        function p = getProperty(self, pname)
            nProps = length(self.m_listGuiProperties);
            % build strings
            for i=1:nProps
                p = self.m_listGuiProperties(i);
                if strcmp(p.getName(), pname)
                    return
                end
            end
        end      
        %------------------------------------------------------------------
        function self = resetGuiProps(self)
            childs = get(self.m_hguiPropertyGrid, 'Children');
            for i=1:length(childs)
                delete(childs(i));
            end
        end 
        %------------------------------------------------------------------
        function CBLoadButton(self, hObject, eventdata)     
            [FileName,PathName] = uigetfile('*.props.mat','Select the config file'); 
            if ~isempty(FileName) && ~isnumeric(FileName)
                D = load(fullfile(PathName, FileName));
                for i=1:length(D.props)
                    p = self.getProperty(D.props{i,1});
                    p.set(D.props{i,2});
                end
            end
        end        
        %------------------------------------------------------------------
        function CBSaveButton(self, hObject, eventdata)     
            [FileName,PathName] = uiputfile('*.props.mat','Select the config file');  
            if ~isempty(FileName) && ~isnumeric(FileName)
                props = self.getPropertyCell();
                date_ = date();
                version = 1;
                readme = 'This file was created with the Matlab property tools. Dont edit if you dont know what you are doing.';
                save(fullfile(PathName, FileName), 'props', 'date_', 'version', 'readme', '-v7.3');
            end  
        end
    end
end