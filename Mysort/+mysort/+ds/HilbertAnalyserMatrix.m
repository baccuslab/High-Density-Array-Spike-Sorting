classdef HilbertAnalyserMatrix < mysort.ds.FrequencyAnalyserMatrix
    properties

    end
    methods
        %------------------------------------------------------------------
        function self = HilbertAnalyserMatrix(X, samplesPerSecond, hpf, lpf, order, varargin)
            self = self@mysort.ds.FrequencyAnalyserMatrix(X, samplesPerSecond, hpf, lpf, order, varargin);
         
        end
        %------------------------------------------------------------------
        function X = getData_(self, idx1, idx2)
            X = getData_@mysort.ds.FrequencyAnalyserMatrix(self,idx1, idx2);
            X = sum(abs(hilbert(X)),2);
        end
        %------------------------------------------------------------------
        function varargout = size(self, varargin)
            dims = [size(self.X,1) 1];
            varargout = matlabfilecentral.parseSize.parseSize(dims, nargout, varargin{:}); 
        end           
    end
end