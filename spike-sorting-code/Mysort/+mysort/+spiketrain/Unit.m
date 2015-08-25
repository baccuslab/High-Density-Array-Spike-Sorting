classdef Unit < handle
    
    properties
        parentSpikeSorting
        spikeTrains
        ID
        footprint
        amplitude
    end
    
    methods
        % -----------------------------------------------------------------
        function self = Unit(ID, st, parentSpikeSorting, footprint)
            if nargin < 4
                footprint = [];
                if nargin < 3
                    parentSpikeSorting = [];
                    if nargin < 2
                        st = [];
                        if nargin < 1
                            ID = [];
                        end
                    end
                end
            end
                  
            self.parentSpikeSorting = parentSpikeSorting;
            self.ID = ID;
            self.spikeTrains = st;
            self.footprint = footprint;
            self.parentSpikeSorting = parentSpikeSorting;
            self.amplitude = max(abs(footprint(:)));
        end
        
        % -----------------------------------------------------------------
        function plot(self)
            
        end
        
        % -----------------------------------------------------------------
        function st = getSpikeTrains(self)
            st = self.spikeTrains;  
        end        
    end
end