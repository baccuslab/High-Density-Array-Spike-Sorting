classdef ChipView < handle
    % To rotate 180
    % set(gca, 'XDir', 'reverse', 'YDir', 'normal')
    properties  
        % Gui elements
        Window
        ChipAxes
        ChipAxesBG
        callback
        
        % internal vars
        config
        chipSpecification
        bPlotInitialized
        axisLimsMinMax
        Properties
        highdensEventdata
        backgroundImage
        
        % rectangle vars
        startPoint
        endPoint
        rectangle
        
        % plotting Handles
        m_hElectrodes
        m_hSelectedElectrodes
        m_hTextLabels
        m_hRoutedElectrodes
        m_hBackground
    end
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = ChipView(windowHandle, axesHandle, axesHandleBG, callback, chipSpecification, props)
            self.Window = windowHandle;
            self.ChipAxes = axesHandle;
            self.ChipAxesBG = axesHandleBG;
            self.callback = callback;
            self.config = [];
            self.Properties = props;
            self.setChipSpecification(chipSpecification);
        end
        %------------------------------------------------------------------
        function setChipSpecification(self, spec)
            self.chipSpecification = spec;

            self.bPlotInitialized = false;
            self.config = [];
            self.reset();
        end
        %------------------------------------------------------------------
        function setConfig(self, config)
            self.config = config;
            self.plot();
        end
        %------------------------------------------------------------------
        function setSingleChannelDetectedSpikes(self, TCA, CL)
            display('not implemented');
%             self.update();
        end
        %------------------------------------------------------------------
        function reset(self)
            self.rectangleReset();
            self.startPoint = [];
            self.endPoint = [];
            self.rectangle = [];
            self.m_hElectrodes = [];
            self.m_hSelectedElectrodes = [];
            self.m_hTextLabels = [];
            self.m_hRoutedElectrodes = [];  
%             set(self.ChipAxes, 'xtick', [], 'xticklabel', [xsp]);
%             set(self.ChipAxes, 'ytick', [], 'yticklabel', []);            
            if ~isempty(self.chipSpecification)
                xmima = [min(self.chipSpecification.x) max(self.chipSpecification.x)];
                ymima = [min(self.chipSpecification.y) max(self.chipSpecification.y)];
                xslack = .1*diff(xmima);
                yslack = .1*diff(ymima);
                self.axisLimsMinMax = [xmima+[-1 1]*xslack ymima+[-1 1]*yslack];
            else
                self.axisLimsMinMax = [0 1 0 1];
            end
            set([self.ChipAxes self.ChipAxesBG], 'xlim', self.axisLimsMinMax(1:2), ...
                               'ylim', self.axisLimsMinMax(3:4));    
            self.plot();
        end
        %------------------------------------------------------------------
        function setHighDensEventDataBackground(self, event_data)
            self.highdensEventdata = event_data;
%             self.backgroundImage = [];
            self.plot();
        end
        %------------------------------------------------------------------
        function setBackgroundImage(self, img)
%             self.highdensEventdata = [];
            self.backgroundImage = img;
            self.plot();
        end        
        %------------------------------------------------------------------
        function resetBackground(self)
            self.backgroundImage = [];
            self.highdensEventdata = [];
            self.plot();
        end
        %------------------------------------------------------------------
        function plot(self)
            bothAx = [self.ChipAxes self.ChipAxesBG];
            self.axisLimsMinMax(1:2) = get(self.ChipAxes, 'xlim');
            self.axisLimsMinMax(3:4) = get(self.ChipAxes, 'ylim');
            
%             set(self.ChipAxes, 'HitTest', 'off');
            cla(self.ChipAxes);
            if isempty(self.chipSpecification)
                return
            end
            % make axis ij
            set(bothAx,'XDir','normal','YDir','reverse');
            % make axis equal
            axis(bothAx, 'equal');    
            set(bothAx, 'nextplot', 'add');
            bShowBackground = self.Properties.getPropertyValue('showBackground');
            bZoomOnBackground = false;
            if bShowBackground && ~(isempty(self.highdensEventdata) && isempty(self.backgroundImage))
                bZoomOnBackground = self.Properties.getPropertyValue('autoZoomOnBackground');
                self.plotBackground(bZoomOnBackground);
                backgroundColorMap = self.Properties.getPropertyValue('colorMap');
                colormap(backgroundColorMap);
            end            

            self.plotElectrodes();
            
%             for i=1:size(ELPOS,1)
%                 text(ELPOS(i,1), ELPOS(i,2), num2str(ELPOS(i,3)), 'parent', self.ChipAxes);
%             end

            self.plotSelectedElectrodes();
            bPlotLabels = self.Properties.getPropertyValue('showTextLabels');
            if bPlotLabels
                self.plotTextLabels();
            end
            self.plotRoutedElectrodes();
            self.plotStimElectrodes();
            if ~bZoomOnBackground
%                 set(self.ChipAxes, 'xtick', [], 'xticklabel', []);
%                 set(self.ChipAxes, 'ytick', [], 'yticklabel', []);
                set(bothAx, 'xlim', self.axisLimsMinMax(1:2), ...
                                   'ylim', self.axisLimsMinMax(3:4));            
                self.bPlotInitialized = true;
            end
            bRotate180 = self.Properties.getPropertyValue('rotate180');
            % axis ij,             'normal' 'reverse'
            % axis xy,             'normal' 'normal'
            % axis ij + rotate180, 'reverse' 'normal'
            if bRotate180
                set(bothAx, ...
                    'XDir','reverse',...
                    'YDir','normal');
            else
                set(bothAx, ...
                    'XDir','normal',...
                    'YDir','reverse');
            end
            xlabel(self.ChipAxes, 'X [\mum]');
            ylabel(self.ChipAxes, 'Y [\mum]');
            box(self.ChipAxes, 'on');
            box(self.ChipAxesBG, 'off');
            linkaxes(bothAx, 'xy');

%             axis(self.ChipAxes, 'tight');
            
            
%             if strcmp(self.lastAction, 'add')
%                 selEl = ELPOS(self.currentElectrodes, 1:2);
%                 plot(self.ChipAxes, selEl(:,1), selEl(:,2), 'om', 'markersize', 20, 'linewidth', 2);
%             end
        end
        %------------------------------------------------------------------
        function plotBackground(self, bZoomOnBackground)
            if isempty(self.highdensEventdata) && isempty(self.backgroundImage)
                return
            end
            bothAx = [self.ChipAxes self.ChipAxesBG];
            if nargin ==1
                bZoomOnBackground = false;
            end
            axes(self.ChipAxes);
            if nr.isEventMap(self.highdensEventdata)
                hidens_generate_event_map(self.highdensEventdata, 'no_new_fig',...
                    'freqthr', 0.1, 'markerplot', 'border', 15,...
                    'markerplotfreq', 'dist_scale', ...
                    'dist_scale_len', 200, 'scalebar_simple');            
            elseif nr.isAlignedImage(self.highdensEventdata)
                % This is Davids aligned microscope images
                imagesc([self.highdensEventdata.xstart self.highdensEventdata.xend], ...
                        [self.highdensEventdata.ystart self.highdensEventdata.yend], ...
                        self.highdensEventdata.IM);
                if bZoomOnBackground
                    set(bothAx, 'xlim', [.98*self.highdensEventdata.xstart 1.02*self.highdensEventdata.xend],...
                                       'ylim', [.98*self.highdensEventdata.ystart 1.02*self.highdensEventdata.yend]);
                end
            end
            if ~isempty(self.backgroundImage)
                if isstruct(self.backgroundImage)
                else
                    imagesc(self.backgroundImage);
                end
            end
%             self.m_hBackground = plot(self.ChipAxes, ...
%                 1000, ...
%                 1000, 'kx', 'markersize', 30);
        end
        
        %------------------------------------------------------------------
        function plotElectrodes(self)
            marker = self.Properties.getPropertyValue('electrodeMarker');
            markerSize = self.Properties.getPropertyValue('electrodeMarkerSize');
            notDummyEl = self.chipSpecification.dummy==0;
            self.m_hElectrodes = plot(self.ChipAxes, ...
                self.chipSpecification.x(notDummyEl), ...
                self.chipSpecification.y(notDummyEl), marker, 'markersize', markerSize);          
        end
        %------------------------------------------------------------------
        function plotSelectedElectrodes(self)
            if isempty(self.config) || isempty(self.config.selectedElectrodes)
                return
            end
            
            selectedEl = self.config.selectedElectrodes(:,1); 
            selectedPriority = self.config.selectedElectrodes(:,2); 
            marker = self.Properties.getPropertyValue('selectionMarker');
            priorities  = [50 -1 0];
            markerSizes = [5 7 9];
            self.m_hSelectedElectrodes =  [];
            for i=1:length(priorities)
                idx = find(selectedPriority == priorities(i));
                tmp = plot(self.ChipAxes, ...
                self.chipSpecification.x(selectedEl(idx)), ...
                self.chipSpecification.y(selectedEl(idx)), marker, 'markersize', markerSizes(i), 'linewidth', 2);
                self.m_hSelectedElectrodes = [self.m_hSelectedElectrodes; tmp];
            end
        end
        %------------------------------------------------------------------        
        function plotTextLabels(self)
            if isempty(self.config) || isempty(self.config.selectedElectrodes)
                return
            end
            selectedEl = self.config.selectedElectrodes(:,1);                      
            self.m_hTextLabels = zeros(1, size(selectedEl,1));
            for i=1:length(selectedEl)
                self.m_hTextLabels(i) = text(...
                        self.chipSpecification.x(selectedEl(i)) +2, ...
                        self.chipSpecification.y(selectedEl(i)),...
                        num2str(self.chipSpecification.el_idx(selectedEl(i))), 'parent', self.ChipAxes);
            end
        end
        %------------------------------------------------------------------        
        function plotRoutedElectrodes(self)
            if ~isfield(self.config, 'routedIdx')
                return
            end       
            selectedEl = self.config.routedIdx;
            marker = self.Properties.getPropertyValue('routedMarker');
            bPlotChannelText = self.Properties.getPropertyValue('plotChannelText');
            offsetx = 8;
            offsety = 8;
            if bPlotChannelText && isfield(self.config,'routedReadOutChannels')
                for i=1:length(selectedEl)
                    self.m_hTextLabels(end+1) = text(self.chipSpecification.x(selectedEl(i))+offsetx, ...
                             self.chipSpecification.y(selectedEl(i))+offsety,...
                             num2str(self.config.routedReadOutChannels(i)), 'parent', self.ChipAxes);
                end
            end
            self.m_hRoutedElectrodes = plot(self.ChipAxes, ...
                self.chipSpecification.x(selectedEl), ...
                self.chipSpecification.y(selectedEl), marker, 'markersize', 8, 'linewidth', 2);            
        end        
        %------------------------------------------------------------------        
        function plotStimElectrodes(self)
            if ~isfield(self.config, 'stimElIdx')
                return
            end       
            selectedEl = self.config.stimElIdx;
            marker = self.Properties.getPropertyValue('stimElMarker');
            self.m_hRoutedElectrodes = plot(self.ChipAxes, ...
                self.chipSpecification.x(selectedEl), ...
                self.chipSpecification.y(selectedEl), marker, 'markersize', 16, 'linewidth', 1);            
        end            
        %------------------------------------------------------------------
        function self = CBButtonDown(self, hObject, eventdata)
            [point relPoint] = nr.getClickPoint(self.ChipAxes);
            self.rectangleReset();
            if any(relPoint<0) || any(relPoint>1) 
                return
            end
            self.startPoint = point(1,1:2);
        end        
        %------------------------------------------------------------------
        function self = CBMouseMove(self, hObject, eventdata)       
            if isempty(self.startPoint)
                return
            end
            point = nr.getClickPoint(self.ChipAxes);
            self.endPoint = point;
            if ishandle(self.rectangle)
                delete(self.rectangle);
            end
            x1 = self.startPoint(1);  
            x2 = point(1,1);
            y1 = self.startPoint(2);
            y2 =  point(1,2);
            x = min(x1, x2);
            y = min(y1, y2);
            w = max(.0001,abs(x1-x2));
            h = max(.0001,abs(y1-y2));
            position = [x y w h];
            
            self.rectangle = rectangle('Position', position, 'parent', self.ChipAxes);
        end
        %------------------------------------------------------------------
        function self = CBButtonUp(self, hObject, eventdata)
            if isempty(self.startPoint)
                self.rectangleReset();
                return
            end
            if isempty(self.endPoint)
                self.endPoint = self.startPoint + [1 1];
            end
            
            p1 = self.startPoint;
            p2 = self.endPoint;
            x1 = min(p1(1), p2(1));
            x2 = max(p1(1), p2(1));
            y1 = min(p1(2), p2(2));
            y2 = max(p1(2), p2(2));
            self.rectangleReset();
            self.callback(x1,x2,y1,y2);
        end

        %------------------------------------------------------------------
        function rectangleReset(self)
            self.startPoint = [];
            self.endPoint = [];
            if ishandle(self.rectangle)
                delete(self.rectangle);
            end
            self.rectangle = [];
        end
    end
end