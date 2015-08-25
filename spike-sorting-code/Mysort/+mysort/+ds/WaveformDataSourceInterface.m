classdef WaveformDataSourceInterface < handle
    properties
        name
        MultiElectrode
        samplesPerSecond
        SpikeSortingContainers
        activeSpikeSortingIdx
        bReturnSortingResiduals
    end
    
    methods(Abstract)
        getWaveform(self, t, cutLeft, cutLength, channelidx)
    end
    
    methods
        %------------------------------------------------------------------
        function self = WaveformDataSourceInterface(name, samplesPerSecond, MultiElectrode, SpikeSortingContainers)
            if nargin > 2 && ~isempty(MultiElectrode)
                assert(isa(MultiElectrode, 'mysort.ds.MultiElectrode'), 'MultiElectrode must be a mysort.ds.MultiElectrode!');
                self.MultiElectrode = MultiElectrode;
            else
                self.MultiElectrode = [];
            end
            assert(isempty(samplesPerSecond) || isnumeric(samplesPerSecond), 'invalid samplerate!');
            assert(length(samplesPerSecond)<2, 'samplerate must not be a vector!');
            
            self.bReturnSortingResiduals = 0;
            self.samplesPerSecond = double(samplesPerSecond);
            self.name = name;
            if nargin == 4
                self.SpikeSortingContainers = SpikeSortingContainers;
                self.activeSpikeSortingIdx = 1;
            else
                self.SpikeSortingContainers = {};
                self.activeSpikeSortingIdx = [];
            end
        end
        %------------------------------------------------------------------        
        function sr = getSamplesPerSecond(self)
            sr = self.samplesPerSecond;
        end
        %------------------------------------------------------------------        
        function sr = getSampleRate(self)
            sr = self.getSamplesPerSecond();
        end
        %------------------------------------------------------------------
        function si = getSampleInterval(self)
            if isempty(self.samplesPerSecond)
                si = [];
                return
            end
            si = 1/self.samplesPerSecond;
        end 
        %------------------------------------------------------------------
        function b = hasSpikeSorting(self)
            b = ~isempty(self.SpikeSortingContainers);
        end
        %------------------------------------------------------------------
        function S = getActiveSpikeSorting(self)
            S = [];
            idx = self.activeSpikeSortingIdx;
            if isempty(idx)
                if self.hasSpikeSorting()
                    self.activeSpikeSortingIdx = 1;
                    idx = 1;
                else
                    return
                end
            end
            S = self.getSpikeSorting(idx);
        end
        %------------------------------------------------------------------
        function S = getActiveSpikeSortingContainer(self)
            idx = self.activeSpikeSortingIdx;
            S = [];
            if isempty(idx)
                return
            end            
            S = self.getSpikeSortingContainer(idx);
        end
        
        %------------------------------------------------------------------
        function found = getSpikeSortingIdx4Name(self, name)
            found = [];
            for i=1:length(self.SpikeSortingContainers)
                if ~isempty(self.SpikeSortingContainers{i}) && ...
                    strcmp(self.SpikeSortingContainers{i}.name, name)
                    found = i;
                    return
                end
            end
            warning('could not find the requested name!');
        end
        %------------------------------------------------------------------
        function idx = setActiveSpikeSorting(self, name)
            idx = self.getSpikeSortingIdx4Name(name);
            assert(~isempty(idx), 'No spike sorting with that name!');
            self.activeSpikeSortingIdx = idx;
        end
        
        %------------------------------------------------------------------
        function setActiveSpikeSortingIdx(self, idx)
            if ischar(idx) 
                self.setActiveSpikeSorting(idx);
            else
                assert(idx>0 && idx <= length(self.SpikeSortingContainers), 'Index out of bounds!');
                self.activeSpikeSortingIdx = idx;
            end
        end
        
        %------------------------------------------------------------------
        function S = getSpikeSorting(self, idx)
            if isempty(idx) || isempty(self.SpikeSortingContainers)
                S = [];
                return
            end
            S = self.SpikeSortingContainers{idx};
        end
        %------------------------------------------------------------------
        function S = getSpikeSortingContainer(self, idx)
            if isempty(idx) || isempty(self.SpikeSortingContainers)
                S = [];
                return
            end
            S = self.SpikeSortingContainers{idx};
        end        
        %------------------------------------------------------------------
        function addSpikeSorting(self, S)
            self.SpikeSortingContainers{end+1} = S;
            if isempty(self.activeSpikeSortingIdx)
                self.activeSpikeSortingIdx = 1;
            end
        end
    end
end