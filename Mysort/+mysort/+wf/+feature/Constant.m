classdef Constant < mysort.wf.feature.BufferedInterface
    properties             
       
    end
    
    methods
        %------------------------------------------------------------------
        function self = MinMax(wfm)
            self = self@mysort.wf.feature.BufferedInterface('MinMax', wfm);
        end
        
        %------------------------------------------------------------------
        function F = calc_(self, wfs)
            F = [min(wfs(:, 3:15), [], 1) max(wfs(:, 3:15), [], 1)];
        end        
    end
end