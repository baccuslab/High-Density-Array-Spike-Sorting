classdef ElectrodeFeature < hdmeagui.Feature
    properties
        ME
        elnumber
    end
    methods
        %------------------------------------------------------------------
        function self = ElectrodeFeature(name, eventTimes, eventIDs, eventFeatures, ME, elnumber)
            assert(isa(ME, 'mysort.ds.MultiElectrode'), 'ME must be a MultiElectrode!');
            self = self@hdmeagui.Feature(name, eventTimes, eventIDs, eventFeatures);
            self.ME = ME;
            self.elnumber = elnumber;
        end
    end
end