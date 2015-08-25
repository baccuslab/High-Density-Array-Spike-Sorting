
classdef Chunker < handle
 properties      
        P
        bHasNextChunk
        currentChunk
        nChunks
        Len
        chunks
        
        chunk_start
        chunk_end
        chunk_start_overlap
        chunk_end_overlap
        
        progress
    end
    
    methods (Abstract)
    end
    
    methods
        %%% ------------------------------------------------------
        % TODO: MIN CHUNK SIZE. 
        function self = Chunker(Len, varargin)
            self.P.chunkSize = 5000000;
            self.P.overlap = 0; % obsolete, use chunkOverlap
            self.P.chunkOverlap = 0;
            self.P.offset = 0;
            self.P.minChunkSize = 1;
            self.P.progressDisplay = 'none'; % or 'console' or 'progressbar';
            self.P.showTotalTimeAtEnd = 1;
            self.P = mysort.util.parseInputs(self.P, varargin, 'error');            
            assert(isa(Len, 'double'), 'Len must be of type double, otherwise weird things can happen, e.g., with Matlabs min function!');
            assert(isa(self.P.chunkSize, 'double'), 'self.P.chunkSize must be of type double, otherwise weird things can happen, e.g., with Matlabs min function!');
            assert(isa(self.P.chunkOverlap, 'double'), 'isa(self.P.chunkOverlap must be of type self.P.chunkOverlap, otherwise weird things can happen, e.g., with Matlabs min function!');
            assert(isa(self.P.offset, 'double'), 'self.P.offset must be of type self.P.offset, otherwise weird things can happen, e.g., with Matlabs min function!');
            assert(isa(self.P.minChunkSize, 'double'), 'self.P.minChunkSize must be of type double, otherwise weird things can happen, e.g., with Matlabs min function!');
            assert(isnumeric(Len), 'Len must be numeric!');
            assert(Len > 0, 'Len must be >0!');
            assert(Len >= self.P.minChunkSize, 'minChunkSize is bigger than Len!');
            
            if ~self.P.overlap == 0
                warning('In Chunker the variable overlap was set. This is obsolete. Use chunkOverlap instead!')
                self.P.chunkOverlap = self.P.overlap;
            end
            
            self.Len = Len;
            self.progress = [];
            self.reset();
        end
        
        %%% ------------------------------------------------------
        function reset(self)
            self.precomputeChunks();
            %self.Len = self.P.stopOffset - (self.P.startOffset-1); 
            self.nChunks = size(self.chunks,1);
            self.currentChunk = 1;
            self.bHasNextChunk = true;
            
            self.chunk_start   = [];
            self.chunk_end     = [];
            self.chunk_start_overlap = [];
            self.chunk_end_overlap   = []; 
            
            if strcmp(self.P.progressDisplay, 'none')
                self.progress = [];
            elseif strcmp(self.P.progressDisplay, 'console')
                self.progress = mysort.util.ProcessTimer(self.nChunks);
            elseif strcmp(self.P.progressDisplay, 'progressbar')
                matlabfilecentral.progressbar.progressbar;
            else
                error('Unknown progress display (%s)', self.P.progessDisplay);
            end            
        end
               
        %%% ------------------------------------------------------
        function [chunkOverlap chunk chunkLen] = getNextChunk(self)
            assert(self.hasNextChunk()==true, 'No more chunks available!');
          
            self.chunk_start = self.chunks(self.currentChunk,1);
            self.chunk_end   = self.chunks(self.currentChunk,2);
            
            self.chunk_start_overlap = max(1, self.chunk_start-self.P.chunkOverlap);
            self.chunk_end_overlap   = min(self.Len, self.chunk_end +self.P.chunkOverlap);            
            self.currentChunk = self.currentChunk +1;
            if self.chunk_end >= self.Len
                self.bHasNextChunk = false;
            end
            chunkOverlap = [self.chunk_start_overlap self.chunk_end_overlap];
            chunk = [self.chunk_start self.chunk_end];
            chunkLen = chunk(2)-chunk(1)+1;
            
            chunkOverlap = chunkOverlap + self.P.offset;
            chunk = chunk + self.P.offset;
            
            if strcmp(self.P.progressDisplay, 'none')
                self.progress = [];
            elseif strcmp(self.P.progressDisplay, 'console')
                self.progress.next(self.currentChunk-1);
                self.progress.showProgress();
                if self.P.showTotalTimeAtEnd && ~self.hasNextChunk()
                    self.progress.showTotal();
                end
            elseif strcmp(self.P.progressDisplay, 'progressbar')
                matlabfilecentral.progressbar.progressbar((self.currentChunk-1)/self.nChunks);
            else
                error('Unknown progress display (%s)', self.P.progessDisplay);
            end        

        end
        
        %%% ------------------------------------------------------
        function b = hasNextChunk(self)
            b = self.bHasNextChunk;
        end
        %%% ------------------------------------------------------
        function precomputeChunks(self)
            L = self.Len;
            cs = self.P.chunkSize;
            mincs = self.P.minChunkSize;
            
            % Check if we have only one chunk
            if L <= cs
                self.chunks = [1 L];
                return
            end
            
            self.chunks = [(0:cs:(L-cs))+1; cs:cs:L]';
            % Check if Length is multiple of chunksize
            if round(L/cs) == L/cs
                return
            end
                
            
            % try to make last chunk shorter than all others
            self.chunks(end+1,:) = [self.chunks(end,2)+1 L];
            % check if last chunk is still long enough
            lastChunkLen = self.chunks(end,2)-self.chunks(end,1)+1;
            if lastChunkLen >= mincs
                return
            end
            
            % last chunk is too short. lets try to make the last two chunks
            % shorter but still longer than the minchunksize
            missingChunkLen = mincs - lastChunkLen;
            secondLastChunkLenAfterCorrection = cs - missingChunkLen;
            if secondLastChunkLenAfterCorrection > mincs
                % is long enough, so correct
                self.chunks(end-1, 2) = self.chunks(end-1,2)-missingChunkLen;
                self.chunks(end,   1) = self.chunks(end-1, 2)+1;
                return
            end
            error('The minimal chunksize does not allow to split last two chunks. Increase chunkSize or decrease minChunkSize');
        end
    end
end