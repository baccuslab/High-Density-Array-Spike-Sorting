
classdef DataHandle < mysort.datafile.DataFileInterface
    properties
        X
    end
    
    methods
        %%% ------------------------------------------------------
        function self = DataHandle(X, samplesPerSecond, varargin)
            if nargin == 1
                samplesPerSecond = 20000;
            end
            self = self@mysort.datafile.DataFileInterface(X, samplesPerSecond, varargin{:}); %'DataHandle',
        end
        %%% ------------------------------------------------------
        function [bLoaded nC Len] = init(self, X, samplesPerSecond, varargin)
%             assert(isnumeric(X), 'X must be numerical!');
            assert(~isempty(X), 'X must not be empty!');
            assert(ndims(X) == 2, 'X must have 2 dimensions! rows = recording channels');
            assert(samplesPerSecond == round(samplesPerSecond), 'samplesPerSecond must be an integer!');
            if nargin < 3
                samplesPerSecond = 20000;
            end
            self.X = X;
            bLoaded = true;
            nC  = size(X,1);
            Len = size(X,2);
            self.samplesPerSecond = samplesPerSecond;
        end
        
        %%% ------------------------------------------------------
        function X = getData_(self, start, stopp, channels)
            X = self.X(channels, start:stopp);
        end 
    end
end