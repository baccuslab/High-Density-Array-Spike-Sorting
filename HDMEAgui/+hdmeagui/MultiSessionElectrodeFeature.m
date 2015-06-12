classdef MultiSessionElectrodeFeature < hdmeagui.ElectrodeFeature
    properties
        eventSessions
    end
    methods
        %------------------------------------------------------------------
        function self = MultiSessionElectrodeFeature(name, eventTimes, eventIDs, eventSessions, eventFeatures, ME, elnumber)
            assert(isa(ME, 'mysort.ds.MultiSessionMultiElectrode'), 'ME must be a MultiSessionMultiElectrode!');
            self = self@hdmeagui.ElectrodeFeature(name, eventTimes, eventIDs, eventFeatures, ME, elnumber);
            self.eventSessions = eventSessions;
        end
        %------------------------------------------------------------------
        function c = getElectrodeCardinality(self)
            c = self.ME.getElectrodeCardinality(self.elnumber);
        end
        %------------------------------------------------------------------
        function removeCloseEvents(self, blockEventsMSGdf, slack)
            for i=1:length(blockEventsMSGdf)
                gdf = blockEventsMSGdf{i};
                if isempty(gdf)
                    continue
                end
                sidx = find(self.eventSessions == i);
                if isempty(sidx)
                    continue;
                end
                M = mysort.spiketrain.findNearSpikes(self.eventTimes(sidx), gdf(:,2), slack);
                ridx = sidx(M(:,1));
                self.eventSessions(ridx) = [];
                self.eventTimes(ridx) = [];
                self.eventIDs(ridx) = [];
                self.eventFeatures(ridx) = [];
            end
        end
    end
end