classdef SortingGuiElement < guiutil.GuiElement
    properties
        % gui elements
        DataSource
        
        TemplateListBox
        SortingListBox
        selectedSortingContainer
%         currentWaveformManager
        
        hfuncSortingSelectionCB
        hfuncTemplateSelectionCB
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = SortingGuiElement(varargin)
            P.sortingCallback = [];
            P.templateCallback = [];
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);
            self = self@guiutil.GuiElement(uP{:});   
            self.makeLayout();
            self.hfuncSortingSelectionCB = P.sortingCallback;
            self.hfuncTemplateSelectionCB = P.templateCallback;
        end
        %------------------------------------------------------------------
        function makeLayout(self)
            layout = uiextras.VBoxFlex( 'Parent', self.Parent, 'Spacing', 3 );

            self.SortingListBox = uicontrol( 'Style', 'list', ...
               'Max', 1,...
               'BackgroundColor', 'w', ...
               'Parent', layout, ...
               'String', '', ...
               'Selected', 'off',...
               'Value', 1, ...
               'Callback', @(x,y) self.CBSortingListBoxSelection(x,y));
            template_box = uiextras.BoxPanel( ...
               'Parent', layout, ...
               'Title', 'Templates');
                self.TemplateListBox = uicontrol( 'Style', 'list', ...
                   'Max', 1000,...
                   'BackgroundColor', 'w', ...
                   'Parent', template_box, ...
                   'String', '', ...
                   'Selected', 'off',...
                   'Value', 1, ...
                   'Callback', @(x,y) self.CBTemplateListBoxSelection(x,y));
            set(layout, 'Sizes', [-1 -4]);   
        end      
           
        %------------------------------------------------------------------
        function self = setDataSource(self, ds)
            self.DataSource = ds;
            self.update();
        end
        
        %------------------------------------------------------------------
        function MEA = getMEA(self)
            MEA = self.MEA;
        end        

        %------------------------------------------------------------------
        function self = update(self)
            if isempty(self.DataSource)
                self.reset();
                return
            end
            v_old = get(self.SortingListBox, 'Value');
            SC = self.DataSource.SpikeSortingContainers;
            sort_str = {};
            for i=1:length(SC)
                if isempty(SC{i}.name)
                    sort_str{i} = 'no_name';
                else
                    sort_str{i} = SC{i}.name;
                end
            end
            if isempty(v_old) | v_old > length(sort_str)
                v = 1;
            else
                v = v_old;
            end
            set(self.SortingListBox, 'String', sort_str, 'Value', v, 'Selected', 'off');
%             if v ~= v_old;
                self.setSortingIndex(v);
%             end
        end
        %------------------------------------------------------------------
        function reset(self)
            set(self.TemplateListBox, 'String', '', 'Value', 1);
            set(self.SortingListBox, 'String', '', 'Value', 1);
        end
        %------------------------------------------------------------------
        function setSortingIndex(self, idx)
            if length(self.DataSource.SpikeSortingContainers)<idx
                set(self.TemplateListBox, 'String', '');
                return
            end
            SC = self.DataSource.SpikeSortingContainers{idx};
            self.selectedSortingContainer = SC;
            t_str = {};
            for i=1:length(SC.unitNames)
                t_str{i} = num2str(SC.unitNames(i));
            end
            set(self.TemplateListBox, 'String', t_str);
            if get(self.TemplateListBox, 'Value') > length(t_str)
                set(self.TemplateListBox, 'Value', length(t_str));
            end
            if isempty(get(self.TemplateListBox, 'Value')) || get(self.TemplateListBox, 'Value') <= 0 
                set(self.TemplateListBox, 'Value', 1);
            end
        end
        
        %------------------------------------------------------------------
        function self = CBSortingListBoxSelection(self, hObject, eventdata) 
            v = get(self.SortingListBox, 'Value');
            if ~isempty(v)
                self.setSortingIndex(v)
            end
            if ~isempty(self.hfuncSortingSelectionCB)
                self.hfuncSortingSelectionCB(v);
            end
        end
        %------------------------------------------------------------------
        function self = CBTemplateListBoxSelection(self, hObject, eventdata) 
            v = get(self.TemplateListBox, 'Value');
            if ~isempty(self.hfuncTemplateSelectionCB)
                self.hfuncTemplateSelectionCB(v);
            end
        end        
    end
end
