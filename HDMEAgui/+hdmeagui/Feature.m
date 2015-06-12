classdef Feature < handle
    properties
        name
        eventTimes
        eventIDs
        eventFeatures
    end
    methods
        %------------------------------------------------------------------
        function self = Feature(name, eventTimes, eventIDs, eventFeatures)
            self.name = name;
            assert(length(eventTimes) == length(eventIDs), 'an event must have a time, and id and a feature!');
            assert(length(eventTimes) == length(eventFeatures), 'an event must have a time, and id and a feature!');
            self.eventTimes = eventTimes;
            self.eventIDs = eventIDs;
            self.eventFeatures = eventFeatures;
        end
        %------------------------------------------------------------------
        function idx = getIdx4Threshold(self, thr)
            if thr>=0
                idx = find(self.eventFeatures > thr);
            else
                idx = find(self.eventFeatures < thr);
            end
        end     
        %------------------------------------------------------------------
        function id = getIDs4Threshold(self, thr)
            idx = self.getIdx4Threshold(thr);
            id = self.eventIDs(idx);
        end  
    end
end