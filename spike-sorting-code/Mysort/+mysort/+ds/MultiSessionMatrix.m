classdef MultiSessionMatrix < mysort.ds.MultiSessionInterface
    properties

    end
    
   
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = MultiSessionMatrix(name, X_cell, s_per_sec)
            if ~iscell(X_cell)
                X_cell = {X_cell};
            end
            nSessions = length(X_cell);
            SL = mysort.ds.Matrix.empty();
            for i=1:nSessions
                SL(i) = mysort.ds.Matrix(X_cell{i}, s_per_sec);
            end
            self = self@mysort.ds.MultiSessionInterface(name, s_per_sec, SL);
        end
    end
end