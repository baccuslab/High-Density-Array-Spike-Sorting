classdef TemplateInterface < mysort.wf.feature.Interface
    properties
       template
    end
    
    %------------------------------------------------------------------
    methods (Abstract)
        get(self)
    end
    
    methods
        %------------------------------------------------------------------
        function self = TemplateInterface(name, wfm, template)
            self = self@mysort.wf.feature.Interface(name, wfm);
            self.template = template;
        end
    end
end