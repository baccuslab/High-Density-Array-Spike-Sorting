classdef SliderAxes < guiutil.SingleAxes
    properties
        % gui elements
        mainLayout
        Slider
        
        SliderPositionEditProperties
        SetSliderPositionButton
        
        WindowEditProperties
        SetWindowWitdhButton
        
        % gui internals
        sliderPosition
        sliderPositionName
        sliderPositionDisplayName
        bAvoidReentrant
        bAvoidKeyReentrant
        
        windowWidthMax
        windowWidthMin        
    end
    
 
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = SliderAxes(varargin)
            P.name = 'position';
            P.displayName = 'Position';
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);       
            self = self@guiutil.SingleAxes(uP{:});
            self.bAvoidReentrant = false;
            self.bAvoidKeyReentrant = false;
            self.sliderPosition = 1;
            self.sliderPositionName = P.name;
            self.sliderPositionDisplayName = P.displayName;
        end
        %------------------------------------------------------------------
        function makeLayout(self, nAxes)
            if nargin == 1
                nAxes = 1;
            end
            p = self.Parent;
            try
                set(p, 'KeyPressFcn', @self.CBKeyPressed);
            catch
            end
            ml = uiextras.VBox( 'Parent', p, 'Spacing', 3, 'BackgroundColor', 'w');
                for i=1:nAxes     
                    self.Axes(i) = axes('Parent', ml);
                           
                    hhAxes = handle(self.Axes(i));  % hAxes is the Matlab handle of our axes
                    hProp = findprop(hhAxes,'xlim');  % a schema.prop object
                    hListener = handle.listener(hhAxes, hProp, 'PropertyPostSet', @self.CBAxesXLimChange);
                    setappdata(self.Axes(i), 'XLimListener', hListener);
                end
                sliderLayout = uiextras.HBox( 'Parent', ml, 'Spacing', 3 );
                    self.Slider = uicontrol('Parent', sliderLayout,...
                        'Style','slider','Min',0,'Max',1,...
                        'Value',.5, 'SliderStep', [0 1],...
                        'KeyPressFcn', @self.CBKeyPressed,...
                        'callback', []); %@self.CBSliderChange
                    hhSlider = handle(self.Slider);
                    hProp = findprop(hhSlider, 'Value');  
                    hListener = handle.listener(hhSlider,hProp,'PropertyPostSet',@self.CBSliderPositionChange);
                    setappdata(self.Slider, 'ValueListener', hListener);
                    
                    setTimeButtonLayout = uiextras.HBox( 'Parent', sliderLayout, 'Spacing', 3);
                        props = {self.sliderPositionName, 1, 'displayName', self.sliderPositionDisplayName};
                        self.SliderPositionEditProperties = guiutil.GuiProperties(setTimeButtonLayout, props,...
                            'loadSaveButton', 0, 'stringLayoutWidth', 40); 
%                         self.SetSliderPositionButton = uicontrol( 'Style', 'PushButton', ...
%                            'Parent', setTimeButtonLayout, ...
%                            'String', 'Go', ...
%                            'Callback', @self.CBSetSliderPositionButton);   
                        props = {'window', 1000, 'displayName', 'window'};
                        self.WindowEditProperties = guiutil.GuiProperties(setTimeButtonLayout, props,...
                            'loadSaveButton', 0, 'stringLayoutWidth', 40); %, 'callback', @self.CBSliderWindowChange
                        self.SetWindowWitdhButton = uicontrol( 'Style', 'PushButton', ...
                           'Parent', setTimeButtonLayout, ...
                           'String', 'Go', ...
                           'Callback', @self.CBSetWindowWidthButton);   

                    set(setTimeButtonLayout, 'Sizes', [115, 115, 40]);
                set(sliderLayout, 'Sizes', [-1, 276]);
            set(ml, 'Sizes', [-1*ones(1, nAxes) 20]);
            self.mainLayout = ml;
            self.bAvoidReentrant = false;
        end
        
        %------------------------------------------------------------------
        function setSliderPosition(self, pos, updateSlider)
            if self.bAvoidReentrant 
                return
            end
%             caller_name = mysort.util.getCallerName();
%             disp(['setSliderPosition - ' caller_name])  
            self.bAvoidReentrant = true;
            self.sliderPosition = pos;
            p = self.SliderPositionEditProperties.getProperty(self.sliderPositionName);
            p.set(pos);
            if updateSlider
                set(self.Slider, 'Value', pos);
            end
%             self.setAllListeners('off')
            self.update();
%             self.setAllListeners('on');            
            self.bAvoidReentrant = false;
        end
        %------------------------------------------------------------------
        function setAllListeners(self, state)
            % seems not to work. why?
            hListener = getappdata(self.Slider, 'ValueListener');
            set(hListener, 'Enabled', state);
            for i=1:length(self.Axes)
                hListener = getappdata(self.Axes(i), 'XLimListener');
                set(hListener, 'Enabled', state);
            end
        end
        %------------------------------------------------------------------
        function setXLim(self, xlim)
            if self.bAvoidReentrant 
                return
            end            
            w = self.WindowEditProperties.getProperty('window');
            w.set(xlim(2)-xlim(1));
            self.setSliderPosition(xlim(1),1);
        end        
%         %------------------------------------------------------------------
%         function CBSliderChange(self, a, b)
%             if self.bAvoidReentrant 
%                 return
%             end            
%             t = get(self.Slider, 'Value');
%             self.setSliderPosition(t,0);            
%         end        
        %------------------------------------------------------------------
        function CBAxesXLimChange(self, hProp, eventData)
            if self.bAvoidReentrant 
                return
            end            
           hAxes = eventData.AffectedObject;
           xlim = get(hAxes,'XLim');
           self.setXLim(xlim);
        end 
        %------------------------------------------------------------------
        function CBSetSliderPositionButton(self, a, b)
            self.bAvoidReentrant = false;
            p = self.SliderPositionEditProperties.getProperty('time');
            t = p.get();
            self.setSliderPosition(t,1);
        end        
        %------------------------------------------------------------------
        function CBSliderPositionChange(self, hProp, eventData)
            if self.bAvoidReentrant 
                return
            end            
%             caller_name = mysort.util.getCallerName();
%             disp(['CBSliderPositionChange - ' caller_name])           
            hAxes = eventData.AffectedObject;
            t = get(hAxes,'Value');
            self.setSliderPosition(t,0);             
        end 
        %------------------------------------------------------------------
        function CBSetWindowWidthButton(self, hProp, eventData)
            self.bAvoidReentrant = false;
            p = self.SliderPositionEditProperties.getProperty('time');
            t = p.get();
            self.setSliderPosition(t, 1);    
        end
        %------------------------------------------------------------------
        function CBKeyPressed(self, a, evt)
            if self.bAvoidKeyReentrant
                return
            end
            self.bAvoidReentrant = false;
            self.bAvoidKeyReentrant = true;
            w = self.WindowEditProperties.getProperty('window');
            p = self.SliderPositionEditProperties.getProperty('time');
            sl_min = get(self.Slider, 'Min');
            sl_max = get(self.Slider, 'Max');
            sl_maxW = sl_max-sl_min;
            switch evt.Key
                case 'uparrow'
                    oldWidth = w.get();
                    center = p.get() + oldWidth/2;
                    newWidth = max(self.windowWidthMin, min(self.windowWidthMax, min(oldWidth*2, sl_maxW)));
                    newStart = center - newWidth/2;
                    w.set(newWidth);
                    p.set(newStart);
                case 'downarrow'
                    oldWidth = w.get();
                    center = p.get()+oldWidth/2;
                    newWidth = max(self.windowWidthMin, min(self.windowWidthMax, min(oldWidth/2, sl_maxW)));
                    newStart = center - newWidth/2;
                    w.set(newWidth);
                    p.set(newStart);
                case 'leftarrow'
                    t = p.get();
                    p.set(max(sl_min, t-w.get()/2));
                case 'rightarrow'
                    t = p.get();
                    p.set(min(sl_max, t+w.get()/2));
                case 'home'
                    p.set(sl_min);
                case 'end'
                    p.set(sl_max-w.get());
                case 'pageup'
                case 'pagedown'
            end            
            self.setSliderPosition(p.get(), 1); 
            self.bAvoidKeyReentrant = false;
        end     
        %------------------------------------------------------------------
%         function CBSliderWindowChange(self, hProp, eventData)
%             self.bAvoidReentrant = false;
%             w = self.WindowEditProperties.getProperty('window');
% %             p = self.SliderPositionEditProperties.getProperty('time');
%             sl_min = get(self.Slider, 'Min');
%             sl_max = get(self.Slider, 'Max'); 
%             sl_maxW = sl_max-sl_min;
%             set(self.Slider, 'SliderStep', [w.get()/2 w.get()]/sl_maxW);
%         end        
       
    end
end