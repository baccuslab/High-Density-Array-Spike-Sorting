
classdef ComplexDataInterface < mysort.util.DataFileInterface
    properties             
        path
        multiElectrodes
        multiElectrodeDescription
        files
        fileDescription
        interFileDistance
    end
    
    methods (Abstract)
        [bLoaded nC Len] = init(self, fname);
        X = getData_(self, start, stopp);
    end
    
    methods
        %%% ------------------------------------------------------
        function self = ComplexDataInterface(fname, varargin)
            self = self@mysort.util.DataFileInterface(varargin{:});
            self.P.chunkSize = 5000000;
            self.P = mysort.util.parseInputs(self.P, 'DataFileInterface', varargin);            
            self.fname = fname;
            [self.bLoaded self.nC self.Len] = self.init(fname);
          
        end
        
        %%% ------------------------------------------------------
        function X = getData(self, varargin)
            assert(self.bLoaded==true, 'File was not sucessfully initialized!');
            p.start = [];
            p.stopp = [];
            p = mysort.util.parseInputs(p, 'DataFileInterface.getData', varargin);
            if isempty(p.start); p.start=1; end 
            if isempty(p.stopp); p.stopp=self.Len; end
            
            assert(p.start > 0, 'start must be greater than 0');
            assert(p.stopp <= self.Len, 'stopp must be smaller than Len!'); 
            
            X = self.getData_(self.P.startOffset-1+p.start, self.P.startOffset-1 + p.stopp);
        end        
    end
    
    methods (Static)
    end    
end