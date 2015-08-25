classdef binaryFileMatrix < handle
    
    % Creates a matlab style matrix from a binary file and/or write data
    % into this file.
    % A 2D matrix is stored column after column in the file:
    %           | a b c |
    % M(:,:) =  | d e f |      --> M(:) = [a d b e c f]
    
    properties (Constant)
        
    end
    properties
        fname
        memmap
        dims
        P
    end
    
    methods
        %% CONSTRUCTOR FUNCTION
        %------------------------------------------------------------------
        function self = binaryFileMatrix(fname, dims, varargin)
            self.fname = fname;
            self.dims = double(dims);
            
            P.writable = false;
            P.precision  = 'int16';
            P.offset = 0;
            
            self.P = mysort.util.parseInputs(P, varargin, 'error');
            
            if exist(fname, 'file')
                self.openMemmap();
            else
                assert(self.P.writable, 'File does not exist, cannot open readonly!');
                self.initializeFile();
                self.openMemmap();
            end
            
        end
        
        %------------------------------------------------------------------
        function openMemmap(self)
            self.memmap = memmapfile(self.fname, 'Format', {self.P.precision self.dims([2 1]) 'X'}, 'Writable', self.P.writable, 'Offset', self.P.offset);
            assert( size(self.memmap.Data,1)==1, 'Given dimensions don''t match the size of the binary file!')
            %             if  size(self.memmap.Data,1) ~= 1
            %                 warning('Given dimensions don''t match the size of the binary file!');
            %                 self.dims(2) = self.dims(2)*size(self.memmap.Data,1);
            %                 self.memmap = memmapfile(self.fname, 'Format', {self.P.precision self.dims 'X'}, 'Writable', self.P.writable, 'Offset', self.P.offset);
            %             end
        end
        
        %------------------------------------------------------------------
        function initializeFile(self)
            fileID = fopen(self.fname,'a+'); fwrite(fileID, zeros(self.dims), self.P.precision); fclose(fileID);
        end
        
        % DESTRUCTOR
        %------------------------------------------------------------------
        %function delete(self)
        %
        %end
        
        %% TOP LEVEL DATA ACCESS FUNCTIONS
        %------------------------------------------------------------------
        function B = subsref(self,S)
            if strcmp(S(1).type, '.')
                % for the '.' subsref call the standard one
                B = builtin('subsref', self, S);
                return
            end
            assert(strcmp(S.type, '()'), 'Only () is implemented!');
            assert(length(S)==1, 'Only X(a,b) is possible, no further subindexing!');
            assert(length(S.subs) < 4, 'Only 3 dimensions can be indexed!');
            
            B = self.getData(S.subs{:});
        end
        
        %------------------------------------------------------------------
        function self = subsasgn(self, S, B)
            if strcmp(S(1).type, '.')
                % for the '.' subsref call the standard one
                B = builtin('subsasgn', self, S, B);
                return
            end
            assert(strcmp(S.type, '()'), 'Only () is implemented!');
            assert(length(S)==1, 'Only X(a,b) is possible, no further subindexing!');
            assert(length(S.subs) < 3, 'Only 2 dimensions can be indexed!');
            
            self.setData(B,S.subs{:});
        end
        
        %------------------------------------------------------------------
        function L = getNSamples_(self)
            L = self.dims(1);
        end
        
        %------------------------------------------------------------------
        function varargout = size(self,varargin)
            varargout = matlabfilecentral.parseSize.parseSize( self.dims ,nargout,varargin{:});
        end
        
        %------------------------------------------------------------------
        function out = end(self, k, n)
            out = self.dims(k);
            if n == 1
                out = prod(self.dims);
            end
        end
        %------------------------------------------------------------------
        function self = transpose(self)
            error('not implemented, conflicts with handle class');
        end
        %------------------------------------------------------------------
        function self = ctranspose(self)
            error('not implemented, conflicts with handle class');
        end
        
        %% MIDDLE LEVEL ACCESS FUNCTIONS
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function X = getData(self, varargin)
            if length(varargin) == 2  
                c = {varargin{[2 1]}};
            else
                c = varargin;
            end
            X =  self.memmap.Data.X(c{:})';
            %X =  self.memmap.Data.X(varargin{:})';
        end
        
        %------------------------------------------------------------------
        function setData(self, X, varargin)
            assert( self.P.writable, 'Cannot write to read only file!');
            
            % Make sure that only data with the same number of rows is set:
            dimsX = size(X);
            assert( dimsX(2) == self.dims(2), 'Dimensions missmatch: number of columns must be the same!');
            
            if varargin{1} == ':'
                % Case when call is M(:, :)
                m = min(self.dims(1), dimsX(1));
                v = 1:m;
            else
                % Case when call is M(1:6, :)
                vx = varargin{1};
                assert( length(vx) == dimsX(1), 'Dimensions missmatch: number of rows must be the same!');
                m = min(self.dims(1), vx(end) );
                v = vx(1):m;
            end
            
            self.memmap.Data.X(:,v) = X(1:length(v), :)';
            self.appendData( X((length(v)+1):end, :) );
        end
        
        %------------------------------------------------------------------
        function appendData(self, X)
            if size(X,1) ~= 0
                fileID = fopen(self.fname,'a+');
                fwrite(fileID, X', self.P.precision);
                fclose(fileID);
                self.dims(1) = self.dims(1) + size(X,1);
                self.openMemmap;
            end
        end
        
        %------------------------------------------------------------------
        function wf = getWaveform_(self, nCut, channelindex, cutLength, t1, t2)
%             wf = zeros(nCut, length(channelindex)*cutLength);
            % Build complete index set and access data with a single
            % entry
            IDX = zeros(1, nCut*cutLength);
            for i = 1:nCut
                idx_idx = (i-1)*cutLength +1 : i*cutLength;
                IDX(idx_idx) = t1(i):t2(i); 
            end 
            wf_ = self.getData(IDX, channelindex);
            wf = mysort.wf.t2v(permute(reshape(wf_, [cutLength nCut length(channelindex)]), [1 3 2]));
        end
        
        %------------------------------------------------------------------
        function self = clearFile(self, dims)
            try
                if nargin == 2
                    self.dims = dims;
                end
                self.deleteFile;
                self.initializeFile;
                self.openMemmap;
            catch
                warning(['Could not clear binary file ' self.fname '!']);
            end
        end
        
        %------------------------------------------------------------------
        function self = deleteFile(self)
            try
                delete(self.fname);
            catch
                warning(['Could not delete binary file ' self.fname '!']);
            end
        end
        
    end
end