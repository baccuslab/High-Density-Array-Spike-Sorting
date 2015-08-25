
classdef DataSource < mysort.datasource.DataSourceInterface
    properties
        X
    end
    
    methods
        %------------------------------------------------------------------
        function self = DataSource(X, varargin)
            % Input defs
            self = self@mysort.datasource.DataSourceInterface(varargin{:});
            
            assert(~isempty(X), 'handled data cant be empty!');
%             assert(isnumeric(X), 'X must be numerical!');
            assert(~isempty(X), 'X must not be empty!');
            assert(ndims(X) == 2, 'X must have 2 dimensions! rows = recording channels');
            
            % Construct
            self.X = X;
        end
        %------------------------------------------------------------------        
        function X = getData(self, epochs, channels)
            % Return X but only in the epochs given in "epochs". epochs is
            % a matrix with two columns. First column indicates the start
            % sample of an epoch, the second the end sample.
            % If your data source provides a clever way to implement this
            % overwrite this method in your inherited class implementation
            if ~exist('epochs', 'var') || isempty(epochs)
                epochs = [1 self.getLen()];
            end
            if ~exist('channels', 'var') || isempty(channels)
                channels = 1:self.getNChannel();
            end            
            if size(epochs,1) == 1
                 % this is faster for the easy case
                 X = self.X(channels, epochs(1,1):epochs(1,2));
            else
                epochLengths = mysort.epoch.length(epochs);
                totalEpochLength = sum(epochLengths);
                X = zeros(length(channels), totalEpochLength);
                epochStartIdx = 1;
                for i=1:size(epochs,1)
                    X(:, epochStartIdx:epochStartIdx+epochLengths(i)-1) = ...
                        self.X(channels, epochs(i,1):epochs(i,2));
                    epochStartIdx = epochStartIdx + epochLengths(i);
                end
            end
        end
        %------------------------------------------------------------------        
        function Len = getLen(self)
            Len = size(self.X,2);
        end
        %------------------------------------------------------------------        
        function nC = getNChannel(self) 
            nC = size(self.X, 1);
        end
    end
end