
classdef DataFileInterface < mysort.util.DebuggableClass
    properties             
        fname
        bLoaded
        bHasNextChunk
        currentChunk
        nChunks
        nC
        Len
        samplesPerSecond
    end
    
    methods (Abstract)
        [bLoaded nC Len] = init(self, fname, varargin);
        X = getData_(self, start, stopp, channels);
    end
    
    methods
        %%% ------------------------------------------------------
        function self = DataFileInterface(fname, samplesPerSecond, varargin)
            self = self@mysort.util.DebuggableClass(varargin{:});
            self.P.chunkSize = 5000000;
            self.P.startOffset = [];
            self.P.stoppOffset = [];
            self.P.verbose = self.LEVEL_CHUNKING;
            [self.P remainingParas] = mysort.util.parseInputs(self.P, 'DataFileInterface', varargin, 1);            
            self.fname = fname;
            [self.bLoaded self.nC self.Len] = self.init(fname, samplesPerSecond, mysort.util.deflateP(remainingParas));
            if ~isempty(self.P.startOffset)
            	assert(self.P.startOffset>0, 'startOffset must lie within the datafile!');
                assert(self.P.startOffset<=self.Len, 'startOffset must lie within the datafile!');
            else
                self.P.startOffset = 1;
            end
            if ~isempty(self.P.stoppOffset)
            	assert(self.P.stoppOffset>0, 'stoppOffset must lie within the datafile!');
                assert(self.P.stoppOffset<=self.Len, sprintf('stoppOffset must lie within the datafile! (Len = %d, offset=%d)', self.Len, self.P.stoppOffset));
            else
                self.P.stoppOffset = self.Len;
            end            
            assert(self.P.startOffset < self.P.stoppOffset, 'startOffset must be smaller than stoppOffset');
            
            self.Len = self.P.stoppOffset - (self.P.startOffset-1); 
            self.nChunks = ceil(self.Len/self.P.chunkSize);
            if self.bLoaded
                self.bHasNextChunk = true;
                self.currentChunk = 1;
            end
        end
        
        %%% ------------------------------------------------------
        function [X start stopp] = getData(self, varargin)
            assert(self.bLoaded==true, 'File was not sucessfully initialized!');
            p.start = [];
            p.stopp = [];
            p.channels = [];
            p = mysort.util.parseInputs(p, 'DataFileInterface.getData', varargin);
            if isempty(p.start); p.start=1; end 
            if isempty(p.stopp); p.stopp=self.Len; end
            
            assert(p.start > 0, 'start must be greater than 0');
            assert(p.start <= self.Len, 'start must be smaller than Len!'); 
            assert(p.stopp > 0, 'stopp must be greater than 0');
            assert(p.stopp <= self.Len, 'stopp must be smaller than Len!'); 
            if isempty(p.channels); p.channels = 1:self.nC; end
            assert(ndims(p.channels)==2, 'channels must be a vector of channel indices!');
            assert(min(p.channels)>=1, 'min channel index must be larger than zero!');
            assert(max(p.channels)<=self.nC, 'max channel index must be smaller or equal to channel number!');
            
            X = self.getData_(self.P.startOffset-1+p.start, self.P.startOffset-1 + p.stopp, p.channels);
            start = p.start;
            stopp = p.stopp;
        end
        
        %%% ------------------------------------------------------
        function [X chunkOffset] = getNextChunk(self)
            assert(self.bLoaded==true, 'File was not sucessfully initialized!');
            assert(self.hasNextChunk()==true, 'No more chunks available!');
            
            chunkStart = self.P.chunkSize*(self.currentChunk-1) +1;
            chunkStopp = min(chunkStart+self.P.chunkSize-1, self.Len);
            channels = 1:self.nC;
            X = self.getData_(self.P.startOffset-1 + chunkStart,  self.P.startOffset-1 + chunkStopp, channels);
            self.currentChunk = self.currentChunk +1;
            if chunkStopp >= self.Len
                self.bHasNextChunk = 0;
            end
            chunkOffset = chunkStart-1;
        end
        
        %%% ------------------------------------------------------
        function b = hasNextChunk(self)
            b = self.bHasNextChunk;
        end
    end
    
    %%% -------------------------------------------------------------------
    methods (Static)
    end    
end