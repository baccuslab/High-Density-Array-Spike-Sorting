classdef Gui < handle
    properties(Constant=true)
        EL_FIXED = 1;
        EL_HIDENS = 2;
        EL_HIGHPRIORITY = 3;
        EL_MOVING = 4;
    end
    properties
        % gui elements
        Window
        mainLayout
        ControlTabPanel
        MEAElement
        ViewPanel
        controlPanel
        ListBox
        CopyButton
        RenameButton
        LoadButton
        SaveButton
        RouteButton
        RouteAllButton
        LoadBackgroundButton
        LoadLastPicButton
        ResetBackgroundButton
        SendConfigButton
        StartRecordingButton
        
        SelModeButton
        PriorityModeDropDown
        SelStimModeButton
        ResetConfigButton
        
        PropertyPanel
        
        ChipAxes
        ChipAxesBackground
        renameEdit
        
        Properties
        ViewGuiProperties
        
        % gui variables
        electrodePriorityCosts
        selectionMode
        selectStimMode
        routeOutputPath
        
        % temp buffers for gui elements
        selectionModeValueBuffer
        
        % Views
        ChipView
        
        % callbacks
        postCallbackCheck
        
        % configs
        chipSpecification
        configFile
        configList
        activeConfig
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function gui = Gui()
            gui.Window = mysort.plot.figure('w', 1000, 'h', 800, 'name', 'NeuroRouter');
            gui.mainLayout = uiextras.HBoxFlex( 'Parent', gui.Window, 'Spacing', 3 );
            gui.ControlTabPanel = uiextras.TabPanel( 'Parent', gui.mainLayout,...
                        'Padding', 5 );
            % Create the panels
            gui.controlPanel = uiextras.BoxPanel( ...
               'Parent', gui.ControlTabPanel, ...
               'Title', 'Select a Config' );
            gui.MEAElement = mysort.mea.H5FileGuiElement(gui.ControlTabPanel,...
                'hasOpenButton', 1,...
                'openCallback', @gui.CBMEAOpenButton,...
                'selectionCallback', @gui.CBMEASessionSelection);
            gui.ControlTabPanel.TabNames = {'Configs', 'MEA File'};
            gui.ControlTabPanel.SelectedChild = 1;      
            
            gui.ViewPanel = uiextras.BoxPanel( ...
               'Parent', gui.mainLayout, ...
               'Title', 'Viewing: ');

            % Adjust the main layout
            set( gui.mainLayout, 'Sizes', [200,-3]  );

            gui.selectionMode = 'add';
            gui.selectStimMode = 'Normal';
            gui.postCallbackCheck = @nr.postCallbackCheck;

            % Create the controls
            % cant put a controlLayout as a member of a struct?!
            controlLayout = uiextras.VBox( 'Parent', gui.controlPanel, ...
               'Padding', 3, 'Spacing', 3 );
            gui.ListBox = uicontrol( 'Style', 'list', ...
               'BackgroundColor', 'w', ...
               'Parent', controlLayout, ...
               'String', '', ...
               'Value', 1, ...
               'Callback', @(x,y) gui.CBListBoxSelection(x,y));

            controlButtonLayout0 = uiextras.HBox( 'Parent', controlLayout, ...
               'Padding', 3, 'Spacing', 3);
            renamestr = uicontrol('Style', 'text', ...
                'Parent', controlButtonLayout0, ...
                'String', 'New Name: ');
            gui.renameEdit = uicontrol('Style', 'edit', ...
                'Parent', controlButtonLayout0, ...
                'String', '');
            
            controlButtonLayout1 = uiextras.HBox( 'Parent', controlLayout, ...
               'Padding', 3, 'Spacing', 3);
            gui.CopyButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', controlButtonLayout1, ...
               'String', 'Copy', ...
               'Callback', @(x,y) gui.CBCopyButton(x,y));
            gui.RenameButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', controlButtonLayout1, ...
               'String', 'Rename', ...
               'Callback', @(x,y) gui.CBRenameButton(x,y));    
           
           
            controlButtonLayout2 = uiextras.HBox( 'Parent', controlLayout, ...
               'Padding', 3, 'Spacing', 3);
            gui.LoadButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', controlButtonLayout2, ...
               'String', 'Load', ...
               'Callback', @(x,y) gui.CBLoadButton(x,y));
            gui.SaveButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', controlButtonLayout2, ...
               'String', 'Save', ...
               'Callback', @(x,y) gui.CBSaveButton(x,y));    

           
            controlButtonLayout3 = uiextras.HBox( 'Parent', controlLayout, ...
               'Padding', 3, 'Spacing', 3);
            gui.RouteButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', controlButtonLayout3, ...
               'String', 'Route', ...
               'Callback', @(x,y) gui.CBRouteButton(x,y));   
            gui.RouteAllButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', controlButtonLayout3, ...
               'String', 'Route All', ...
               'Callback', @(x,y) gui.CBRouteAllButton(x,y));    
           
            controlButtonLayout4 = uiextras.HBox( 'Parent', controlLayout, ...
               'Padding', 3, 'Spacing', 3);
           
           gui.SendConfigButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', controlButtonLayout4, ...
               'String', 'Send Config', ...
               'Callback', @(x,y) gui.CBSendConfigButton(x,y)); 
           
           gui.LoadBackgroundButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', controlButtonLayout4, ...
               'String', 'Load Background', ...
               'Callback', @(x,y) gui.CBLoadBackgroundButton(x,y));
           
            controlButtonLayout5 = uiextras.HBox( 'Parent', controlLayout, ...
               'Padding', 3, 'Spacing', 3);
            
            % load last pic 
            gui.LoadLastPicButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', controlButtonLayout5, ...
               'String', 'Load Last Picture', ...
               'Callback', @(x,y) gui.CBLoadLastPicButton(x,y));
           
            gui.ResetBackgroundButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', controlButtonLayout5, ...
               'String', 'Reset', ...
               'Callback', @(x,y) gui.CBResetBackground(x,y)); 
           
%             gui.StartRecordingButton = uicontrol( 'Style', 'PushButton', ...
%                'Parent', controlButtonLayout5, ...
%                'String', 'Start Record', ...
%                'Callback', @(x,y) gui.CBStartRecordingButton(x,y)); 
%             gui.RouteAllButton = uicontrol( 'Style', 'PushButton', ...
%                'Parent', controlButtonLayout4, ...
%                'String', 'nothing', ...
%                'Callback', @(x,y) x);                
           
            % load property values
            if exist('props.mat', 'file')
                D = load('props.mat');
                props = D.props;
            else
                props = {'showBackground', true, 'displayName', 'Show Background'; ...
                         'colorMap', 'gray', 'displayName', 'BG Colormap'; ...
                         'autoZoomOnBackground', false, 'displayName', 'BG Zoom'; ...
                         'rotate180', true, 'displayName', 'Rotate 180deg'; ...
                         'showTextLabels', false, 'displayName', 'Show Text Labels'; ...
                         'plotChannelText', true, 'displayName', 'Show Routed Labels'; ...
                         'stimElMarker', '.g', 'displayName', 'Stimulation Marker';...
                         'electrodeMarker', '+k', 'displayName', 'Electrode Marker';...
                         'electrodeMarkerSize', 4, 'displayName', 'El. Markersize';...
                         'selectionMarker', 'or', 'displayName', 'Selection Marker'; ...
                         'routedMarker', '.b', 'displayName', 'Routed Marker'};
            end
            cb = @gui.CBPropertyChange;
           
            % create property interface
            gui.PropertyPanel = uiextras.BoxPanel( ...
               'Parent', controlLayout, ...
               'Title', 'Properties');
            gui.Properties = guiutil.GuiProperties(gui.PropertyPanel, props,...
                'callback', cb,...
                'loadSaveButton', 1);           
           
            set(controlLayout, 'Sizes', [-1 28 28 28 28 28 28 280] ); % Make the list fill the space
            
            % Create the view
            props = {'x_spacing', 1, 'displayName', 'x-sp'; ...
                     'y_spacing', 1, 'displayName', 'y-sp'};
            viewLayout = uiextras.VBox( 'Parent', gui.ViewPanel, ...
               'Padding', 2, 'Spacing', 2);
            viewButtonLayout = uiextras.HBox( 'Parent', viewLayout, ...
               'Padding', 2, 'Spacing', 2);
            gui.ViewGuiProperties = guiutil.GuiProperties(viewButtonLayout, props,...
                'callback', [],'loadSaveButton', 0);            
            gui.SelModeButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', viewButtonLayout, ...
               'String', 'Mode: Add', ...
               'Callback', @(x,y) gui.CBSwitchSelModeButton(x, y));     
           
            gui.PriorityModeDropDown = uicontrol( 'Style', 'popup', ...
               'Parent', viewButtonLayout, ...
               'String', 'No Priority|Low Priority|High Priority', ...
               'Callback', @(x,y) gui.CBSwitchPriorityDropdown(x, y));  
           
            gui.SelStimModeButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', viewButtonLayout, ...
               'String', 'Mode: Normal Electrode', ...
               'Callback', @(x,y) gui.CBSwitchStimModeButton(x, y));     
           
            % config reset
            gui.ResetConfigButton = uicontrol( 'Style', 'PushButton', ...
               'Parent', viewButtonLayout, ...
               'String', 'Reset Electrodes', ...
               'Callback', @(x,y) gui.CBResetConfigButton(x, y));   
           
            chipAxesBox = uipanel( 'Parent', viewLayout);
            gui.ChipAxesBackground = axes( 'Parent', chipAxesBox,  'ActivePositionProperty', 'OuterPosition', 'OuterPosition', [0 0 1 1]);
            gui.ChipAxes = axes( 'Parent', chipAxesBox, 'color', 'none', 'ActivePositionProperty', 'OuterPosition', 'OuterPosition', [0 0 1 1]);
            set(viewLayout, 'Sizes', [40 -1] ); % Make the buttons small on top           
            
            gui.chipSpecification = nr.hidens_get_all_electrodes(2);
            gui.ChipView = nr.ChipView(gui.Window, gui.ChipAxes, gui.ChipAxesBackground, @gui.CBacceptRectangle, gui.chipSpecification, gui.Properties);
            gui.configFile = 'configs.mat';
            gui.configList = [];
            gui.routeOutputPath = pwd;            
            
            gui.electrodePriorityCosts = -1;

            % define and register callbacks
            set(gui.Window, 'WindowButtonDownFcn',   @(x,y) gui.ChipView.CBButtonDown(x,y));
            set(gui.Window, 'WindowButtonUpFcn',     @(x,y) gui.ChipView.CBButtonUp(x,y));
            set(gui.Window, 'WindowButtonMotionFcn', @(x,y) gui.ChipView.CBMouseMove(x,y));
            gui.newConfig();
            gui.update();
        end
        
        %------------------------------------------------------------------
        function setActiveConfig(self, cidx)
            assert(cidx <= length(self.configList) && cidx > 0, 'Config index out of bounds!');
            self.activeConfig = cidx;
            self.update();                
        end
        
        %------------------------------------------------------------------
        function newConfig(self)
            c.selectedElectrodes = [];
            c.name = 'unnamedConfig';
            c.routedIdx = [];
            if isempty(self.configList)
                self.configList = c;
            else
                self.configList(end+1) = c;
            end
            self.setActiveConfig(length(self.configList));
        end
        %------------------------------------------------------------------
        function loadConfigFile(self, file)
            if ~exist(file, 'file')
                warning('could not find file: %s', file);
                return
            end
            self.configFile = file;
            D = load(file);
            [pathstr, name, ext] = fileparts(file);
            assert(isfield(D, 'configList'), 'The loaded file was not a valid configList File!');
            assert(isfield(D, 'activeConfig'), 'The loaded file was not a valid configList File!');
            self.configList = D.configList;
            self.activeConfig = D.activeConfig;
            for i=1:length(self.configList)
                self.loadRoutingInformation(pathstr, i);
            end
            self.update();
        end

        %------------------------------------------------------------------
        function update(self)
            cnames = self.getConfigNameStr();
            set(self.ListBox, 'String', cnames, 'value', self.activeConfig);  
            if isempty(self.activeConfig) || isempty(self.configList)
                set(self.ListBox, 'String', '');  
                set(self.renameEdit, 'String', '');
                set(self.ViewPanel, 'Title', 'Viewing: ');
                self.ChipView.reset();
            else
                set(self.ListBox, 'String', cnames, 'value', self.activeConfig);  
                set(self.renameEdit, 'String', self.configList(self.activeConfig).name);
                set(self.ViewPanel, 'Title', sprintf('Viewing: %s', self.configList(self.activeConfig).name));
                self.ChipView.setConfig(self.configList(self.activeConfig));
            end
        end
        %------------------------------------------------------------------
        function str = getConfigNameStr(self)
            str = {};
            for i=1:length(self.configList)
                str{i} = self.configList(i).name;
            end
        end
        
        %------------------------------------------------------------------
        function route(self, pathname, configIdx)  
            self.routeOutputPath = pathname;
            c =  self.configList(configIdx);
            disp('routing');
            disp(c.name);
            if ~isfield(c, 'stimElIdx')
                c.stimElIdx = [];
            end
            npos = nr.hidens_create_route_neuron_cell(c.selectedElectrodes(:,1), c.stimElIdx, c.selectedElectrodes(:,2));
            
            xsp = self.ViewGuiProperties.getPropertyValue('x_spacing');
            ysp = self.ViewGuiProperties.getPropertyValue('y_spacing');
                
            nr.hidens_write_neuroplacement(pathname, [c.name '.neuropos.nrk'], npos,[xsp ysp]);
            nr.hidens_execute_routing(pathname, c.name);
            self.loadRoutingInformation(pathname, configIdx);
            self.update();
        end            
        %------------------------------------------------------------------
        function loadRoutingInformation(self, pathname, configIdx)
            if ~exist(fullfile(pathname, [self.configList(configIdx).name '.el2fi.nrk2']), 'file')
                return
            end
            [routedElNumbers routedReadOutChannels] = nr.hidens_read_el2fi_nrk_file(pathname, self.configList(configIdx).name);
            self.configList(configIdx).routedIdx = self.elNumbers2elIdx(routedElNumbers);
            self.configList(configIdx).routedReadOutChannels = routedReadOutChannels;
        end
        
        %------------------------------------------------------------------
        function idx = elNumbers2elIdx(self, elNumbers)
            idx = zeros(size(elNumbers));
            for i=1:length(elNumbers)
                idx(i) = find(self.chipSpecification.el_idx==elNumbers(i),1);
            end
        end
        
        %------------------------------------------------------------------
        function addElectrodesWithNumbers(self, elnrs)
            idx = self.elNumbers2elIdx(elnrs);
            selectionMode_ = self.selectionMode;
            self.selectionMode = 'add';
            self.acceptElectrodes(idx);
            self.selectionMode = selectionMode_;
        end
        %------------------------------------------------------------------
        function removeElectrodesWithNumbers(self, elnrs)
            idx = self.elNumbers2elIdx(elnrs);
            selectionMode_ = self.selectionMode;
            self.selectionMode = 'remove';
            self.acceptElectrodes(idx);
            self.selectionMode = selectionMode_;
        end        
        
        %------------------------------------------------------------------
        function acceptElectrodes(self, el)
            if isempty(self.activeConfig) || isempty(self.configList)
                self.newConfig();
            end
            config = self.configList(self.activeConfig);     
            if ~isfield(config, 'stimElIdx')
                config.stimElIdx = [];
            end   
            % be backwards compatible, check if the elctrode indices are
            % only a vector, then the priority column is missin
            
            if ~isempty(config.selectedElectrodes) && min(size(config.selectedElectrodes)) == 1 && max(size(config.selectedElectrodes)) > 2
                config.selectedElectrodes = config.selectedElectrodes(:);
                config.selectedElectrodes(:,2) = -1;
            end
            
            dummyIdx = find(self.ChipView.chipSpecification.dummy==1);
            el = setdiff(el, dummyIdx);
            
            if strcmp(self.selectionMode(1:3), 'add')
                if strcmp(self.selectionMode, 'add')
                    % Check which electrodes to select
                    xsp = self.ViewGuiProperties.getPropertyValue('x_spacing');
                    ysp = self.ViewGuiProperties.getPropertyValue('y_spacing');
                    if xsp>1 || ysp > 1
                        elXY = [self.chipSpecification.x(el)' self.chipSpecification.y(el)'];
                        M = nr.makeMatrixFromXY(elXY(:,1), elXY(:,2));
                        % treat 2 rows of M as one row of electrodes. for that we
                        % move every 2 column up by 1.
                        if M(1,1) == 0
                            M = M(1:xsp:end, 2:ysp:end);
                        else
                            M = M(1:xsp:end, 1:ysp:end);
                        end
                        % convert back to indices
                        M_ = M(:);
                        M__ = M_(find(M_));
                        % make subselection
                        el = el(M__);
                    end                
                
                elseif strcmp(self.selectionMode, 'addBlock')
                    % Check which electrodes to select
                    nRows = self.ViewGuiProperties.getPropertyValue('x_spacing');
                    nCols = self.ViewGuiProperties.getPropertyValue('y_spacing');
                    if nRows>1 || nCols > 1
                        colStart = min(self.chipSpecification.x(el));
                        rowStart = min(self.chipSpecification.y(el));
                        elIdx = find(self.chipSpecification.x >= colStart & self.chipSpecification.y >= rowStart);
                        elXY = [self.chipSpecification.x(elIdx)' self.chipSpecification.y(elIdx)'];
                        M = nr.makeMatrixFromXY(elXY(:,1), elXY(:,2));
                        % treat 2 rows of M as one row of electrodes. for that we
                        % move every 2 column up by 1.
                        if M(1,1) == 0
                            s2 = 2;
                        else
                            s2 = 1;
                        end
                        M = M(1:min(end,nRows), s2:min(end,nCols));
                        % convert back to indices
                        M_ = M(:);
                        M__ = M_(find(M_));
                        % make subselection
                        el = elIdx(M__);
                    end                
                else
                    error('Unknown selection mode!')
                end
                % Find electrodes that were already in the list and update
                % priorities
                add = el(:);
                add(:,2) = self.electrodePriorityCosts;
                if ~isempty(config.selectedElectrodes)
                    [C, ia, ib] = intersect(config.selectedElectrodes(:,1), el(:));
                    config.selectedElectrodes(ia,2) = self.electrodePriorityCosts;
                    add(ib,:) = [];
                end
                self.configList(self.activeConfig).selectedElectrodes = [config.selectedElectrodes; add];
                
                if strcmp(self.selectStimMode, 'Stimulation')
                    self.configList(self.activeConfig).stimElIdx = ...
                        union(config.stimElIdx, el);
                end                 
            else 
                % REMOVE electrodes
                [C ia] = setdiff(config.selectedElectrodes(:,1), el);
                self.configList(self.activeConfig).selectedElectrodes = ...
                    self.configList(self.activeConfig).selectedElectrodes(ia,:);
                self.configList(self.activeConfig).stimElIdx = ...
                    setdiff(config.stimElIdx, el);                
            end
            self.ChipView.setConfig(self.configList(self.activeConfig));
        end
        %------------------------------------------------------------------
        function sendConfig(self, fpath, fname)
            if ~isempty(strfind(computer, 'WIN'))
                disp('Sending config only possible under Linux!')
                return
            end
%             fname = 'bloc_exp_002';
            idx = strfind(fname, '.');
            if ~isempty(idx)
                fname = fname(1:idx(1)-1);
            end
            lastpath = pwd;
            cd(fpath);
            unix(sprintf('BitStreamer -n -f %s.hex.nrk2',fname));
            unix(sprintf('nc 11.0.0.7 32124 <%s.cmdraw.nrk2',fname));
            cd(lastpath);
        end        
        %----------------------- MEMBER CALL BACKS ------------------------
        %------------------------------------------------------------------
        function CBacceptRectangle(self, x1, x2, y1, y2)
            el = find(self.chipSpecification.x >= x1 & self.chipSpecification.x <= x2 & ...
                      self.chipSpecification.y >= y1 & self.chipSpecification.y <= y2);
            if isempty(el)
                dx = abs(self.chipSpecification.x - x1);
                dy = abs(self.chipSpecification.y - y1);
                [mi el] = min(dx.^2 + dy.^2);
            end
           
            self.acceptElectrodes(el);
        end
        
        %------------------------------------------------------------------
        function CBPropertyChange(self, p)
            self.update();
        end
        %------------------------------------------------------------------
        function CBMEAOpenButton(self)
            self.CBMEASessionSelection();
        end
        %------------------------------------------------------------------
        function CBMEASessionSelection(self)
            MEA = self.MEAElement.getMEA();
            CL = MEA.getChannelList();
            CL = CL(CL(:,2)==1,3:4)/1000;
%             D = MEA.getPreprocessedData();
%             TCA = D.timesChansAmpsSorted;
%             self.ChipView.setSingleChannelDetectedSpikes(TCA, CL);
        end           
        %----------------------- CALL BACKS -------------------------------
        %------------------------------------------------------------------
        function self = CBListBoxSelection(self, hObject, eventdata)     
            self.setActiveConfig(get(self.ListBox, 'Value'));
        end
        
        %------------------------------------------------------------------
        function CBCopyButton(self, hObject, eventdata)     
            self.configList(end+1) = self.configList(self.activeConfig);  
            self.activeConfig = length(self.configList);
            self.update();
        end        
        %------------------------------------------------------------------
        function CBRenameButton(self, hObject, eventdata)     
            self.configList(self.activeConfig).name = get(self.renameEdit, 'String');
            self.update();
        end        
        
        %------------------------------------------------------------------
        function CBLoadButton(self, hObject, eventdata)     
            [FileName,PathName] = uigetfile('*.nrconfig.mat','Select the config file');   
            self.loadConfigFile(fullfile(PathName, FileName));
        end        
        %------------------------------------------------------------------
        function CBSaveButton(self, hObject, eventdata)     
            [FileName,PathName] = uiputfile('*.nrconfig.mat','Select the config file');   
            % make sure the file ending is correct
            [pathstr, name, ext] = fileparts(FileName);
            FileName = [name '.nrconfig.mat'];
            configList = self.configList;
            activeConfig = self.activeConfig;
            date_ = date();
            version = 1;
            readme = 'This file was created with the Matlab neurorouter. Dont edit if you dont know what you are doing.';
            save(fullfile(PathName, FileName), 'configList', 'activeConfig', 'date_', 'version', 'readme', '-v7.3');
        end   
        %------------------------------------------------------------------
        function CBSwitchSelModeButton(self, hObject, eventdata)
            if isempty(self.selectionModeValueBuffer)
                self.selectionModeValueBuffer = [1 1];
            end
            if strcmp(self.selectionMode, 'add')
                self.selectionMode = 'addBlock';
                set(self.SelModeButton, 'String', 'Mode: Add Block', 'ForegroundColor', [.1 .1 .7]);
                p = self.ViewGuiProperties.getProperty('x_spacing');
                p.setDisplayName('#Columns');
                tmp = p.get();
                p.set(self.selectionModeValueBuffer(1));
                self.selectionModeValueBuffer(1) = tmp;
                
                p = self.ViewGuiProperties.getProperty('y_spacing');
                p.setDisplayName('#Rows');
                tmp = p.get();
                p.set(self.selectionModeValueBuffer(2));
                self.selectionModeValueBuffer(2) = tmp;                
            elseif strcmp(self.selectionMode, 'addBlock')
                self.selectionMode = 'remove';
                set(self.SelModeButton, 'String', 'Mode: Remove', 'ForegroundColor', [.7 .1 .1]);
                p = self.ViewGuiProperties.getProperty('x_spacing');
                p.setDisplayName('n.a.');
                p = self.ViewGuiProperties.getProperty('y_spacing');
                p.setDisplayName('n.a.');   
            else
                self.selectionMode = 'add';
                set(self.SelModeButton, 'String', 'Mode: Add', 'ForegroundColor', [.1 .6 .1]);
                p = self.ViewGuiProperties.getProperty('x_spacing');
                p.setDisplayName('x-sp');
                tmp = p.get();
                p.set(self.selectionModeValueBuffer(1));
                self.selectionModeValueBuffer(1) = tmp;  
                
                p = self.ViewGuiProperties.getProperty('y_spacing');
                p.setDisplayName('y-sp');                 
                tmp = p.get();
                p.set(self.selectionModeValueBuffer(2));
                self.selectionModeValueBuffer(2) = tmp;                  
            end  
        end
        %------------------------------------------------------------------
        function CBSwitchPriorityDropdown(self, hObject, eventdata)
            LOW_PRIORITY = 50;
            NO_PRIORITY = -1;
            HIGH_PRIORITY = 0;
            val = get(hObject, 'Value');
            if val == 1
                self.electrodePriorityCosts = NO_PRIORITY;
            elseif val == 2
                self.electrodePriorityCosts = LOW_PRIORITY; % low
            else
                self.electrodePriorityCosts = HIGH_PRIORITY; % high, 0 is the highest!
            end
        end       
        %------------------------------------------------------------------
        function CBSwitchStimModeButton(self, hObject, eventdata)
            if strcmp(self.selectStimMode, 'Stimulation')
                self.selectStimMode = 'Normal';
                set(self.SelStimModeButton, 'String', 'Electrode Type: Normal');
            else
                self.selectStimMode = 'Stimulation';
                set(self.SelStimModeButton, 'String', 'Electrode Type: Stimulation');
            end  
        end
        %------------------------------------------------------------------
        function CBResetConfigButton(self, hObject, eventdata)
            
            % remove all selected electrodes
            elnrs=0:11015;
            self.removeElectrodesWithNumbers(elnrs)
            
            % self.ChipView.reset() % optionally: reset view
        end
        %------------------------------------------------------------------
        function CBRouteButton(self, hObject, eventdata)     
            [PathName] = uigetdir(self.routeOutputPath, 'Select the config folder for the nrk files');   
            if ~isempty(PathName) && ~isnumeric(PathName)
                self.route(PathName, self.activeConfig);
            end
            self.update();
        end          
        %------------------------------------------------------------------
        function CBSendConfigButton(self, hObject, eventdata)     
            [FileName,PathName] = uigetfile(fullfile(self.routeOutputPath, '*.el2fi.nrk2'), 'Select the nrk file to send');   
            if ~isempty(FileName) && ~isnumeric(FileName)
                self.sendConfig(PathName, FileName);
%                 self.update();                
            end
        end       
%         %------------------------------------------------------------------
%         function CBStartRecordingButton(self, hObject, eventdata)   
%             [PathName] = uigetdir(self.routeOutputPath, 'Select the config folder for the nrk files');   
%             if ~isempty(PathName) && ~isnumeric(PathName)  
%                 self.route(PathName, self.activeConfig);
%                 disp('Start recordings...');
%                 self.update();
%             end
%         end       
        %------------------------------------------------------------------
        function CBRouteAllButton(self, hObject, eventdata)     
            [PathName] = uigetdir(self.routeOutputPath, 'Select the config folder for the nrk files');   
            if ~isempty(PathName) && ~isnumeric(PathName)
                for i=1:length(self.configList)
                    self.route(PathName, i);
                end
                self.update();
            end
        end    
        %------------------------------------------------------------------
        function CBLoadBackgroundButton(self, hObject, eventdata)     
            vars = evalin('base', 'whos()');
            if isempty(vars)
                warning('no variables in global workspace');
                return
            end
            str = {vars.name};
            [s,v] = listdlg('PromptString','Select a high_dense event data variable:',...
                            'SelectionMode','single',...
                            'ListString',str);
            if ~v
                return
            end
            hidens_event_data = evalin('base',vars(s).name);
            if nr.isEventMap(hidens_event_data) || nr.isAlignedImage(hidens_event_data)
                % This is the format the hidense eventmap uses OR
                % This should be Davids format to handle aligned microscope
                % images
                self.ChipView.setHighDensEventDataBackground(hidens_event_data);
            else
                self.ChipView.setBackgroundImage(hidens_event_data);
            end
        end  
                %------------------------------------------------------------------
        function CBLoadLastPicButton(self, hObject, eventdata)     
            
            
            im_path = evalin('base', 'im_path');
            lif_file = evalin('base', 'lif_file');
            cal = evalin('base', 'cal');
            obj_cal = evalin('base', 'obj_cal');
            
                        
            i=get_lif_length([im_path lif_file]);
            dat_last=load_lif_file([im_path lif_file],i);
            pic=transform_align_image(dat_last,1,cal,obj_cal, 'no_plot');

            if nr.isAlignedImage(pic)
                self.ChipView.setHighDensEventDataBackground(pic);
                %self.ChipView.setBackgroundImage(pic);
            end
            
        end  
        %------------------------------------------------------------------
        function CBResetBackground(self, hObject, eventdata)           
            self.ChipView.resetBackground();
        end
    end
end

