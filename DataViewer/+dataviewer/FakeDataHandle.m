classdef FakeDataHandle < dataviewer.DataHandle
    properties
    end
    
    methods
        %------------------------------------------------------------------
        function self = FakeDataHandle(varargin)
            self = self@dataviewer.DataHandle(varargin{:});
            self.P.log_function('This is a fake DB handle with no actual DB connectivity for debug purposes.');
        end
        
        %------------------------------------------------------------------
        function init(self)
            self.P.log_function('This is a fake DB handle with no actual DB connectivity for debug purposes.');
        end
        %------------------------------------------------------------------
        function close(self)
        end
        %------------------------------------------------------------------
        function r = query(self, varargin)
            self.P.log_function('Querying DB...');
            t1 = tic;
            r = {1, 'a', 'b', 'c', 'd'};
            t = toc(t1);
            self.P.log_function(sprintf('... Query done (%.2fs).', t));
        end
    end
end