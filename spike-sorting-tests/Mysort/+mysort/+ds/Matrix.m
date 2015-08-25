classdef Matrix < mysort.ds.ExtendedDataSourceInterface
    properties
        X 
    end
      
    methods
        %------------------------------------------------------------------
        function self = Matrix(X, samplesPerSecond, name, epos, enrs)
            if nargin < 2
                samplesPerSecond = [];
            end
            if nargin < 3
                name = '';
            end
            self = self@mysort.ds.ExtendedDataSourceInterface(name, samplesPerSecond);
            
%            assert(isnumeric(X), 'data must be either a ds.Interface or a matrix!');
            assert(ndims(X)==2, 'data must be either a ds.Interface or a matrix!');  
            nMaxChannels = 2000;
            nMinSamples = 10;
            assert(size(X,2) <= nMaxChannels, sprintf('Too many channels (%d). Do you need to transpose X?', size(X,2)))
            assert(size(X,1) > nMinSamples, sprintf('Too few samples (%d). Do you need to transpose X?', size(X,1)))
            self.X = X;
            if nargin < 5 || isempty(enrs)
                enrs = 1:size(X,2);
            end
            if nargin < 4 || isempty(epos)
                epos = [(1:size(X,2))' zeros(size(X,2),1)];
            end
            self.MultiElectrode = mysort.ds.MultiElectrode(epos, enrs);
        end
        
        %------------------------------------------------------------------
        function X = getData_(self, idx1, idx2)
            X = self.X(idx1, idx2);
        end        
        
        %------------------------------------------------------------------
        function L = getNSamples_(self)
            L = size(self.X,1);
        end
        
        %------------------------------------------------------------------
        function self2 = copy(self)
            self2 = mysort.ds.Matrix(self.X, self.samplesPerSecond, self.name, ...
                self.MultiElectrode.electrodePositions, ...
                self.MultiElectrode.electrodeNumbers);
        end
    end
end