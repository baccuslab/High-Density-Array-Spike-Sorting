classdef Interface < handle
    properties             
       name
       wfm  % waveform manager
    end
    
    %------------------------------------------------------------------
    methods (Abstract)
        get(self)
    end
    
    methods
        %------------------------------------------------------------------
        function self = Interface(name, wfm)
            self.wfm = wfm;
            self.name = name;            
        end
    end
end