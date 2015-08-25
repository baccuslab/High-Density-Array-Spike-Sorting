classdef Projection < mysort.wf.feature.TemplateInterface
    properties
       
    end
    
    methods
        %------------------------------------------------------------------
        function self = Projection(name, wfm, template)
            self = self@mysort.wf.feature.TemplateInterface(name, wfm, template);
        end
        
        %------------------------------------------------------------------
        function P = get(self, idx)
            t = self.template.getWf();
            P = self.wfm.getWfs(idx)*t';
            assert(~any(isnan(P)), 'nan in projections!');
        end        
    end
end