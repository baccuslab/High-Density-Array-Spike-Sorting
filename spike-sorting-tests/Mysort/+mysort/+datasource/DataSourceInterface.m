
classdef DataSourceInterface < mysort.util.DebuggableClass %< handle
    properties
        sampleRate
        sampleInterval        
    end
    
    methods(Abstract)
        getData(self, epochs, channels)
%         getLen(self)
%         getNChannel(self)
    end
    
    methods
        %------------------------------------------------------------------        
        function self = DataSourceInterface(varargin)
            self = self@mysort.util.DebuggableClass(varargin{:});
            self.P.name = [];
            self.P.srate = [];
            self.P.sint = [];
            self.P = mysort.util.parseInputs(self.P, 'DataSource', varargin, 1); 
            
            self.sampleRate = self.P.srate;
            self.sampleInterval = self.P.sint;
            if ~isempty(self.sampleRate) && ~isempty(self.sampleInterval)
                assert((1/self.sampleRate + eps > self.sampleInterval) && ...
                       (1/self.sampleRate - eps < self.sampleInterval), ...
                       'sampleRate is not compatible with sampleInterval');
            end
        end
        
        %------------------------------------------------------------------        
        function sr = getSampleRate(self)
            if ~isempty(self.sampleRate)
                sr = self.sampleRate;
            elseif ~isempty(self.sampleInterval)
                sr = 1/self.sampleInterval;
            else
                sr = 1;
            end
        end
        
        %------------------------------------------------------------------
        function si = getSampleInterval(self)
            if ~isempty(self.sampleInterval)
                si = self.sampleInterval;
            elseif ~isempty(self.sampleRate)
                si = 1/self.sampleRate;
            else
                si = 1;
            end
        end
    end
end