classdef FrequencyAnalyserMatrix < mysort.ds.Matrix
    properties
        filterObject
    end
    methods
        %------------------------------------------------------------------
        function self = FrequencyAnalyserMatrix(X, samplesPerSecond, hpf, lpf, order, varargin)
            self = self@mysort.ds.Matrix(X, samplesPerSecond, varargin{:});
            self.filterObject = mysort.mea.filter_design(hpf, lpf, samplesPerSecond, order);
            assert(ndims(X)==2, 'data must be either a ds.Interface or a matrix!');            
        end
        %------------------------------------------------------------------
        function X = getData_(self, idx1, idx2)
            X = self.X(idx1, idx2);
            X = filtfilthd(self.filterObject, X); 
        end
    end
end