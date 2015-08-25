classdef MCPlot < handle
    properties
        Axes
        dataAxes
        slider
        zoomPopupmenu
        centralMsEdit
        centralMsButton
        zoomWindowCheckbox
        sliderText
    end
    
    methods
        % ------------------------------------------------------
        function self = MCPlot()
            self.Axes = dataViewerAxes();
            children = get(self.Axes, 'Children');
            for i=1:length(children)
                name = get(children(i), 'Tag');
                if strcmp(name, 'dataAxes')
                    self.dataAxes = children(i);
                elseif strcmp(name, 'slider')
                    self.slider = children(i);
                    set(self.slider, 'UserData', @(x) self.sliderMovement(x))
                elseif strcmp(name, 'zoomPopupmenu')
                    self.zoomPopupmenu = children(i);
                    set(self.zoomPopupmenu, 'UserData', @(x) self.zoomPopup(x))
                elseif strcmp(name, 'centralMsEdit')
                    self.centralMsEdit = children(i);
        %             set(self.centralMsEdit, 'UserData', @(x) self.sliderMovement(x))
                elseif strcmp(name, 'centralMsButton')
                    self.centralMsButton = children(i);
                    set(self.centralMsButton, 'UserData', @(x) self.goToMs(x))
                elseif strcmp(name, 'zoomWindowCheckbox')
                    self.zoomWindowCheckbox = children(i);
                    set(self.zoomWindowCheckbox, 'UserData', @(x) self.zoomWindowChecked(x))
                elseif strcmp(name, 'sliderText');
                    self.sliderText = children(i);
                end
            end            
        end
        
        % ------------------------------------------------------
        function bla(self)
            mMax = max(max(m));

            step = 32/numDataPoints * handles.zoomFactor*100;
            handles.sliderStep = step;
            handles.numDataPoints = numDataPoints;

            set(handles.slider,'min',0);
            set(handles.slider,'max',numDataPoints/32);
            set(handles.slider,'sliderStep',[step step]);
            setSliderPosition(handles);
            set(handles.slider,'enable','on');

            set(handles.zoomPopupmenu,'value',1);            
        end
        
        % ------------------------------------------------------
        function zoomWindowChecked(obj, handles)
            obj.updateZoomWindowCheckBox(handles);
            if obj.bZoomWindow
                % draw the zoom window
                obj.drawZoomWindow(handles);    
            else
                % delete the lines that mark the zoom window
                obj.deleteZoomWindow(handles);
            end            
        end
        % ------------------------------------------------------
        function drawZoomWindow(obj,handles)
            temp  = get(handles.dataAxes,'xLim');
            left  = temp(1);
            right = temp(2);
            
            windowBig = right - left;   % the width of the currently viewed area
            windowSmall = windowBig / obj.zoomFactor;   % the width of the area magnified after zooming
            margin = (windowBig - windowSmall) / 2;
            
            left = left + margin;   % the left border of the zoom area
            right = right - margin; % the right border of the zoom area

            % get the upper and lower border of the zoom window
            tmp = get(handles.dataAxes,'ylim');
            lowerB = tmp(1);
            upperB = tmp(2);

            % draw the lines that represent the left and right border
            s = '-.';   % line style
            c = [.8 .8 .8]; % line colour
            h1 = line('xData',[left, left], 'yData',[lowerB, upperB],'Color',c,'LineStyle',s);
            h2 = line('xData',[right, right], 'yData',[lowerB, upperB],'Color',c,'LineStyle',s);

            % store the line handles so the lines can be deleted later
            obj.zoomWindowHandles = [h1 h2];
        end
        % ------------------------------------------------------
        function deleteZoomWindow(obj,handles)
            try
                delete(obj.zoomWindowHandles);
            catch
            end
        end
        % ------------------------------------------------------
        function zoomPopup(obj,handles)
           obj.deleteZoomWindow(handles);

           zf = get(handles.zoomPopupmenu,'value');

           zoomFactor = 1/ (obj.zoomFactor^(zf-1));
           tmp = get(handles.dataAxes,'xLim');
           left   = tmp(1);
           right  = tmp(2);
           window = right-left;
           currentCenter = window/2 + left;
           margin = (obj.numDataPoints/obj.dataHandlerObj.getTetrodeSampleRate(obj.tetrodes(1)) * zoomFactor)/2;
           left = currentCenter - margin;
           right = currentCenter + margin;


           step = obj.sliderStep;
           if step > 0
               step = [step step];
               if (step * zoomFactor) < 1
                   step = step * zoomFactor;
               end
           end
           set(handles.sliderText,'string',['step: ' num2str(step(1)*obj.numDataPoints/obj.dataHandlerObj.getTetrodeSampleRate(obj.tetrodes(1))) 'ms']);
           set(handles.slider,'sliderstep',step);

           set(handles.dataAxes,'xLim',[left right]);
           
           setSliderPosition(handles);
           if get(handles.zoomWindowCheckbox,'value')
               obj.drawZoomWindow(handles);
           end 
        end
        
        % ------------------------------------------------------
        function createSlider(obj,handles)
            [numChannels numDataPoints] = size(obj.data);

            step = 32/numDataPoints * obj.zoomFactor*100;
            obj.sliderStep = step;
            obj.numDataPoints = numDataPoints;

            set(handles.slider,'min',0, 'max',numDataPoints/obj.dataHandlerObj.getTetrodeSampleRate(obj.tetrodes(1)),...
                'sliderStep',[step step]);
            setSliderPosition(handles);
            set(handles.slider,'enable','on');
        end
        
        % ------------------------------------------------------
        function sliderMovement(obj,handles)
            % delete the zoom window first
            obj.deleteZoomWindow(handles);

            pos = get(handles.slider,'Value');
            tmp = get(handles.dataAxes,'xLim');
            left = tmp(1);
            right = tmp(2);
            
            % compute the distance of pos to the edges of the window that is formed by
            % left and right
            leftDistance = pos - left;
            rightDistance = right - pos;
            
            % the slider value should be in the middle of the window, therefore adjust
            % the window (i.e. left and right) accordingly
            meanDist = (leftDistance + rightDistance)/2;
            left = pos - meanDist;
            right = pos + meanDist;

            % set the new window
            set(handles.dataAxes,'xLim',[left right]);
            
            if obj.bZoomWindow
                obj.drawZoomWindow(handles);
            end
        end
        
        % ------------------------------------------------------
        function goToMs(obj,handles)
            obj.updateGoToMs(handles);
            
           % if the user did not provide a valid number string, then end right here
            if isnan(obj.gotoms) 
                return;
            end

            % react if the whised cms is "out of bounds"
            if ((obj.gotoms<0) || (obj.gotoms>obj.numDataPoints))
                return;
            end 
            
            % compute the new xlims around the cms, based on the current xlims
            xlims      = get(handles.dataAxes,'xLim');
            halfWindow = (xlims(2)-xlims(1))/2;
            halfWindow = round(halfWindow);
            xlims(1)   = obj.gotoms - halfWindow;
            xlims(2)   = obj.gotoms + halfWindow;
            
            % set the new xlims 
            set(handles.dataAxes,'xLim',xlims);
            
            % set slider
            tmp   = get(handles.dataAxes,'xLim');
            left  = tmp(1);
            right = tmp(2);
            pos   = (right - left)/2 + left;

            set(handles.slider,'value',pos);
            
            % redraw zoom window
            if obj.bZoomWindow
                obj.deleteZoomWindow(handles);
                obj.drawZoomWindow(handles);
            end
            
        end        
    end
end