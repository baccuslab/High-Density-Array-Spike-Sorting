classdef H5Matrix < mysort.ds.FileBasedMatrix
                    
    properties
        h5DataPath
    end
      
    methods
        %------------------------------------------------------------------
        function self = H5Matrix(fname, h5DataPath, varargin)
            X = mysort.h5.matrix(fname, h5DataPath);
            self = self@mysort.ds.FileBasedMatrix(fname, X, varargin{:});
            self.h5DataPath = h5DataPath;     
        end
    end
end