classdef BufferedInterface < mysort.wf.feature.Interface
    properties             
       fun
       F
       needsUpdate
    end
    
    %------------------------------------------------------------------
    methods (Abstract)
        calc_(self)  % Computes the actual features on the cut spike
                     % waveforms
        getWfs(self) % has to request the spike waveforms from the waveform
                     % manager. In case not all waveforms are needed or
                     % only the waveform of a subset of all channels
    end
    
    methods
        %------------------------------------------------------------------
        function self = BufferedInterface(name, wfm)
            self = self@mysort.wf.feature.Interface(name, wfm);
            self.needsUpdate = true;
        end
        
        %------------------------------------------------------------------
        function calc(self)
            if ~self.needsUpdate
                return
            end
            wfs = self.getWfs();
            self.F = self.calc_(wfs);
            self.needsUpdate = false;            
        end
        %------------------------------------------------------------------
        function F = get(self, idx)
            self.calc();
            F = self.F(idx,:);
        end
    end
end