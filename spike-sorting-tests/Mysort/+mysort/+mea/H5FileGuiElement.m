classdef H5FileGuiElement < handle
    properties
        % gui elements
        Parent
        GlobalDetailPanel
        OpenButton
        fnameStr
        SessionListBox
        DetailPanel
        msgStr
        sizeStr
        fooStr
        ChipElectrodePlot
        
        SortingElement
        
        PreprocessAllButton
        PreprocessButton
        PrefilterAllButton
        PrefilterButton
        
        % callbacks
        hfuncOpenCB
        hfuncSelectCB
        hfuncSortingSelectionCB
        hfuncTemplateSelectionCB
        
        % variables
        hasChipPlot
        hasOpenButton
        hasPreprocessButtons
        fname
        MEA
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = H5FileGuiElement(varargin)
            P.parent = [];
            P.file = [];
            P.hasOpenButton = 1;
            P.hasPreprocessButtons = 1;
            P.hasChipPlot = 1;
            P.openCallback = [];
            P.selectionCallback = [];
            P.sortingSelectionCallback = [];
            P.templateSelectionCallback = [];
            if nargin > 0 && length(varargin{1}) == 1 && (ishandle(varargin{1}) || ...
                              isa(varargin{1}, 'uiextras.Container'))
                P.parent = varargin{1};
                varargin(1) = [];
            end
            P = mysort.util.parseInputs(P, varargin, 'error');
            if nargin==0 || (length(varargin{1}) == 1 && isempty(P.parent) &&  ~(ishandle(varargin{1}) || ...
                              isa(varargin{1}, 'uiextras.Container')))
                P.parent = figure();
            end
            
            self.hasChipPlot = P.hasChipPlot;
            self.hasOpenButton = P.hasOpenButton;
            self.hasPreprocessButtons = P.hasPreprocessButtons;
            self.hfuncOpenCB = P.openCallback;
            self.hfuncSelectCB = P.selectionCallback;
            self.hfuncSortingSelectionCB = P.sortingSelectionCallback;
            self.hfuncTemplateSelectionCB = P.templateSelectionCallback;
            self.Parent = P.parent;            
            self.fname = [];
            self.MEA = [];
            
            self.buildLayout();

            if ~isempty(P.file)
                if ischar(P.file)
                    self.setFileName(P.file);
                elseif isa(P.file, 'mysort.mea.CMOSMEA')
                    self.setMEA(P.file);
                else
                    error('not implemented');
                end
            end
        end
        %------------------------------------------------------------------
        function buildLayout(self)
            mainLayout = uiextras.HBoxFlex( 'Parent', self.Parent, 'Spacing', 3 );
            % LEFT SIDE
            layout = uiextras.VBoxFlex( 'Parent', mainLayout, 'Spacing', 3 );
            sizes = [];
                if self.hasOpenButton
                    self.OpenButton = uicontrol( 'Style', 'PushButton', ...
                       'Parent', layout, ...
                       'String', 'Open', ...
                       'Callback', @self.CBOpenButton);
                   sizes = [sizes 25];
                end

                self.fnameStr = uicontrol( 'Parent', layout, 'Style', 'text', 'String', '');
                sizes = [sizes 30];
                self.SessionListBox = uicontrol( 'Style', 'list', ...
                   'Max', 1000,...
                   'BackgroundColor', 'w', ...
                   'Parent', layout, ...
                   'String', '', ...
                   'Value', 1, ...
                   'Callback', @(x,y) self.CBListBoxSelection(x,y));
               sizes = [sizes -4];
               
                if self.hasPreprocessButtons
                    buttonLayout = uiextras.Grid( 'Parent', layout, 'Spacing', 1 );
                    self.PrefilterButton = uicontrol( 'Style', 'PushButton', ...
                       'Parent', buttonLayout, ...
                       'String', 'Prefilter', ...
                       'Callback', @self.CBPrefilterButton);   
                    self.PrefilterAllButton = uicontrol( 'Style', 'PushButton', ...
                       'Parent', buttonLayout, ...
                       'String', 'Prefilter All', ...
                       'Callback', @self.CBPrefilterAllButton);  
                    self.PreprocessButton = uicontrol( 'Style', 'PushButton', ...
                       'Parent', buttonLayout, ...
                       'String', 'Preprocess', ...
                       'Callback', @self.CBPreprocessButton);  
                    self.PreprocessAllButton = uicontrol( 'Style', 'PushButton', ...
                       'Parent', buttonLayout, ...
                       'String', 'Preprocess All', ...
                       'Callback', @self.CBPreprocessAllButton);                     
                    set(buttonLayout, 'ColumnSizes', [-1 -1], 'RowSizes', [17 17] );
                    sizes = [sizes 35];
                end
                self.DetailPanel = uiextras.BoxPanel( ...
                   'Parent', layout, ...
                   'Title', 'Session Details');
            
                    detailLayout = uiextras.VBox( 'Parent', self.DetailPanel, 'Spacing', 0 );
                        detailGridLayout = uiextras.Grid( 'Parent', detailLayout, 'Spacing', 1 );
                        uicontrol( 'Parent', detailGridLayout, 'Style', 'text', 'String', 'Msg:');
                        uicontrol( 'Parent', detailGridLayout, 'Style', 'text', 'String', 'Size:');
                        uicontrol( 'Parent', detailGridLayout, 'Style', 'text', 'String', 'Foo:');
                        self.msgStr = uicontrol( 'Parent', detailGridLayout, 'Style', 'text', 'String', ':');
                        self.sizeStr = uicontrol( 'Parent', detailGridLayout, 'Style', 'text', 'String', ':');
                        self.fooStr = uicontrol( 'Parent', detailGridLayout, 'Style', 'text', 'String', 'bar');
                        set(detailGridLayout, 'ColumnSizes', [-1 -2], 'RowSizes', [17 17 17] );
                    sizes = [sizes 100];
                    if self.hasChipPlot
                        self.ChipElectrodePlot = mysort.mea.ChipElectrodePlot(detailLayout);
                        set(detailLayout, 'Sizes', [60 -1]);
                        sizes(end) = sizes(end)+200;
                    end            
            set(layout, 'Sizes', sizes);
            % RIGHT SIDE
            sorting_box = uiextras.BoxPanel( ...
                   'Parent', mainLayout, ...
                   'Title', 'Sortings');
            self.SortingElement = mysort.ds.SortingGuiElement('Parent', sorting_box,...
                'sortingCallback', @self.CBSortingListBoxSelection, ...
                'templateCallback', @self.CBTemplateListBoxSelection);
               
            set(mainLayout, 'Sizes', [-7 -3]); 
        end      
           
        %------------------------------------------------------------------
        function self = setFileName(self, fname)
            MEA = mysort.mea.CMOSMEA(fname, 'filterOrder',...
                6, 'hpf', 350, 'lpf', 8000); %2, 'hpf', 500, 'lpf', 3000);
            self.setMEA(MEA);            
        end
        %------------------------------------------------------------------
        function self = setMEA(self, MEA)
            self.MEA = MEA;
            if self.hasChipPlot
                self.ChipElectrodePlot.setMEA(MEA);
            end    
            self.SortingElement.setDataSource(MEA);
            self.update();
        end
        %------------------------------------------------------------------
        function MEA = getMEA(self)
            MEA = self.MEA;
        end        

        %------------------------------------------------------------------
        function self = update(self)
            n = self.MEA.getNSessions();
%             sl = self.MEA.getSessionList();
            s = self.MEA.getActiveSession();
            for i=1:n
                prefilterMarker = '';
                if self.MEA.sessionList(i).isPrefiltered()
                    prefilterMarker = '*';
                end
                preprocessedMarker = '';
                if self.MEA.sessionList(i).isPreprocessed()
                    preprocessedMarker = '+';
                end             
                if self.MEA.sessionList(i).bCovIsPrecalculated
                    preprocessedMarker = [preprocessedMarker 'c'];
                end                        
                sessionstr{i} = [sprintf('Session%d', i) prefilterMarker preprocessedMarker];
            end
            set(self.fnameStr, 'String', self.MEA.fname);
            v = get(self.SessionListBox, 'Value');
            v(v>length(sessionstr)) = [];
            set(self.SessionListBox, 'String', sessionstr, 'Value', v);
            set(self.msgStr, 'String', s.message);
            set(self.sizeStr, 'String', num2str(size(s)));
            if self.hasChipPlot
                self.ChipElectrodePlot.update();
            end            
            if ~isempty(self.SortingElement)
                self.SortingElement.update();
            end
        end

        %------------------------------------------------------------------
        function self = CBListBoxSelection(self, hObject, eventdata)
            if isempty(self.MEA)
                return
            end
            val = get(self.SessionListBox, 'Value');
            if isempty(val)
                return
            end
            self.MEA.setSelectedSessions(val);
            self.update();
            if ~isempty(self.hfuncSelectCB)
                self.hfuncSelectCB();
            end
        end        
        %------------------------------------------------------------------
        function self = CBOpenButton(self, hObject, eventdata)    
            [FileName,PathName] = uigetfile('*.stream.h5','Select the config file');  
            if ~isempty(FileName) & FileName~=0
                self.setFileName(fullfile(PathName, FileName));
                if ~isempty(self.hfuncOpenCB)
                    % if callback is set, call back
                    self.hfuncOpenCB();
                end
            end
        end  
        %------------------------------------------------------------------
        function self = CBPrefilterButton(self, hObject, eventdata)            
            sl = self.MEA.getSelectedSessions();
            for i=1:length(sl)
                if ~sl(i).isPrefiltered()
                    sl(i).prefilter();
                    self.update();
                end
            end
        end
        %------------------------------------------------------------------
        function self = CBPrefilterAllButton(self, hObject, eventdata)            
            sl = self.MEA.getSessionList();
            for i=1:length(sl)
                if ~sl(i).isPrefiltered()
                    sl(i).prefilter();
                    self.update();
                end
            end            
        end
        %------------------------------------------------------------------
        function self = CBPreprocessButton(self, hObject, eventdata)            
            sl = self.MEA.getSelectedSessions();
            for i=1:length(sl)
                if ~sl(i).isPrefiltered()
                    sl(i).prefilter();
                    self.update();
                end
                if ~sl(i).isPreprocessed()
                    sl(i).preprocess();
                    self.update();
                end      
            end              
        end
        %------------------------------------------------------------------
        function self = CBPreprocessAllButton(self, hObject, eventdata)            
            sl = self.MEA.getSessionList();
            for i=1:length(sl)
                if ~sl(i).isPrefiltered()
                    sl(i).prefilter();
                    self.update();
                end
                if ~sl(i).isPreprocessed()
                    sl(i).preprocess();
                    self.update();
                end      
            end              
        end
        %------------------------------------------------------------------
        function self = CBSortingListBoxSelection(self, v) 
            if isempty(self.hfuncSortingSelectionCB)
                return
            end
            if ~isempty(self.MEA)
                self.MEA.setActiveSpikeSortingIdx(v);
            end
            self.hfuncSortingSelectionCB(v);
        end
        %------------------------------------------------------------------
        function self = CBTemplateListBoxSelection(self, v) 
            if isempty(self.hfuncTemplateSelectionCB)
                return
            end  
            self.hfuncTemplateSelectionCB(v);
        end        
    end
end