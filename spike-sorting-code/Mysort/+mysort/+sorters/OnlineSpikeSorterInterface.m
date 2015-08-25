
classdef OnlineSpikeSorterInterface < mysort.sorters.SpikeSorterInterface
    properties       
        chunk_overlap   
        currentChunk
        hasNextChunk
        nChunks
        chunk_start
        chunk_end   
        chunk_start_overlap = [];
        chunk_end_overlap   = [];        
        chunk_sorting
        
        sortBuffer
    end
    
    methods (Abstract)
        adapt(self)
        sortMatrix(self, X)
        plotLastChunkSorting(self)
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------     
        function self = OnlineSpikeSorterInterface(varargin)
            self = self@mysort.sorters.SpikeSorterInterface();
            self.P.chunk_size = 100000;   
            self.P.gui = false;
            self.P = mysort.util.parseInputs(self.P, varargin);
            
            self.chunk_overlap = [];
            self.chunk_start   = [];
            self.chunk_end     = [];
            self.chunk_start_overlap = [];
            self.chunk_end_overlap   = [];
            self.nChunks       = [];
            self.currentChunk  = 0;
            self.hasNextChunk  = 0;
            
            self.chunk_sorting = [];              
        end

        %%% ------------------------------------------------------
%         function b = getHasNextChunk(self)
%             b = self.hasNextChunk;
%         end   
     
        %%% ------------------------------------------------------
        function sorting = sort_(self, X)
            % Cut data into chunks and process them one after the other
            
            self.nChunks = ceil((size(self.DH,1)-self.Tf)/self.P.chunk_size);
            self.chunk_overlap = 2*self.Tf;
            self.hasNextChunk = 1;
            self.currentChunk = 1;
            timer = mysort.util.ProcessTimer(self.nChunks);
            
            sorting = []; 
            while self.hasNextChunk
                self.adapt();
                timer.next();
                str = timer.getProgressString();
                self.debugout(str, self.LEVEL_CHUNKING);
                chunk_sorting = self.sortNextChunk();
                if ~isempty(chunk_sorting)
                    chunk_sorting = sortrows(chunk_sorting,2);
                end
                sorting = [sorting; chunk_sorting];  
                self.sorting = sorting;
                if self.P.gui
                    self.plotLastChunkSorting('gui', 1);
                end
            end
        end
        
        
        %%% ------------------------------------------------------
        function chunk_sorting = sortNextChunk(self) 
            self.chunk_start = (self.currentChunk-1)*self.P.chunk_size + 1;
            self.chunk_end   = min(size(self.DH,1), (self.currentChunk  )*self.P.chunk_size);

            self.chunk_start_overlap = max(1, self.chunk_start-self.chunk_overlap);
            self.chunk_end_overlap   = min(size(self.DH,1), self.chunk_end +self.chunk_overlap);
            
            chunk_sorting = self.sortData(self.chunk_start_overlap, self.chunk_end_overlap);
            
            minRestChunkLen = self.Tf;
            if self.chunk_end_overlap >= size(self.DH,1) - minRestChunkLen
                self.hasNextChunk = 0;
            end   
            
            if ~isempty(chunk_sorting)
                % Shift to current global (not chunk) sample index
                chunk_sorting(:,2) = chunk_sorting(:,2) + self.chunk_start - 1;

                % Clear chunk overlap
                if self.currentChunk > 1
                    % Remove beginning
                    chunk_sorting(:,2) = chunk_sorting(:,2) - (self.chunk_start - self.chunk_start_overlap);
                    chunk_sorting(chunk_sorting(:,2)<self.chunk_start,:) = [];
                end
                if self.currentChunk < self.nChunks
                    % Remove ending
                    chunk_sorting(chunk_sorting(:,2)>self.chunk_end,:) = [];
                end
            end              
            self.chunk_sorting = chunk_sorting;
            self.currentChunk = self.currentChunk + 1;
        end       
        
        %%% ------------------------------------------------------
        function [sorting fOut] = sortData(self, t1, t2) 
            if ~isempty(self.sortBuffer) && self.sortBuffer.t1 == t1 && self.sortBuffer.t2 == t2
                sorting = self.sortBuffer.sorting;
                fOut = self.sortBuffer.fOut;
                return
            end
            CHUNK = self.DH(t1:t2,:)';
                        
            % CALL MAIN SORTING OF AN ONLINE SORTER ON THIS CHUNK
            [sorting fOut] = self.sortMatrix(CHUNK); 
            self.sortBuffer.sorting = sorting;
            self.sortBuffer.fOut = fOut;
            self.sortBuffer.t1 = t1;
            self.sortBuffer.t2 = t2;
        end
    end
end