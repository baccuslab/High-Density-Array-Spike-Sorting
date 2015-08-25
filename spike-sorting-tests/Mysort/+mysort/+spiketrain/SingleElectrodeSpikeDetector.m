 classdef SingleElectrodeSpikeDetector < mysort.spiketrain.SpikeSortingContainerInterface 
    properties
        dataSource
        MultiElectrode
    end
    
    methods(Abstract)
        getUnits(self, unit)
        getSpikeTimes(self, time_idx1, time_idx2, unit)
        size(self, dim)
    end
    
    methods
        %------------------------------------------------------------------
        function self = SingleElectrodeSpikeDetector(dataSource_or_MultiElectrode)
            self = self@mysort.spiketrain.SpikeSortingContainerInterface ('SingleElectrodeSpikeDetector');
            
            if isa(dataSource_or_MultiElectrode, 'mysort.ds.WaveformDataSourceInterface')
                self.dataSource = dataSource_or_MultiElectrode;
                self.MultiElectrode = dataSource_or_MultiElectrode.MultiElectrode;
            else
                self.dataSource = [];
                assert(isa(dataSource_or_MultiElectrode, 'mysort.ds.MultiElectrode'), 'dataSource_or_MultiElectrode must be either of those!')
                self.MultiElectrode = dataSource_or_MultiElectrode;
            end
        end

        %------------------------------------------------------------------
        function getSpikes(self, t1, t2, channels)

        end
        
        %------------------------------------------------------------------
        function [times pks] = detectSpikes(self, varargin)
            [times pks] = mysort.spiketrain.detecSingleElectrodeSpikes(self.dataSource);
        end
    end
 end        