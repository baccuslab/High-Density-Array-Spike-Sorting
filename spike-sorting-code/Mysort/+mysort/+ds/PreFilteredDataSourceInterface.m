classdef PreFilteredDataSourceInterface < mysort.ds.FilteredDataSourceInterface
    properties
        bIsPrefiltered
        prefilterBufferFilename
        prefilterBufferH5info
        prefilterBufferH5path
        prefilterBufferH5matrix
        prefilterScaleFactor
    end
    methods(Abstract)
        getUnfilteredData(self, idx1, idx2)
    end
    
    methods
        %------------------------------------------------------------------
        function self = PreFilteredDataSourceInterface(prefilterBufferFilename, prefilterBufferH5path, h5info, varargin)
            self = self@mysort.ds.FilteredDataSourceInterface(varargin{:});
%             assert(ischar(prefilterBufferFilename), 'prefilterBufferFilename must be a filename!');
%             assert(~isempty(prefilterBufferFilename), 'prefilterBufferFilename must not be empty!');
            self.prefilterBufferFilename = prefilterBufferFilename;
            self.prefilterBufferH5path = prefilterBufferH5path;  
            self.prefilterBufferH5info = h5info;
            self.prefilterScaleFactor = 2^(16-1-1-10)-1; % one for sign, one for range -min:+max, -10 for range from -1024 to 1024 mV
            self.bIsPrefiltered = false;
            
            if exist(prefilterBufferFilename, 'file')
                if isempty(h5info)
                    self.prefilterBufferH5info = hdf5info(prefilterBufferFilename);
                end
                [b h5info Dims bIsArray] = mysort.h5.exist(prefilterBufferFilename, prefilterBufferH5path, self.prefilterBufferH5info);
                % check if prefiltered variable exists and if it is also an
                % array
                if b && bIsArray
                    self.prefilterBufferH5matrix = mysort.h5.matrix(prefilterBufferFilename, prefilterBufferH5path, true);
                    self.bIsPrefiltered = true;
                end
            end
        end
        
        %------------------------------------------------------------------            
        function prefilter(self)
            P.h5chunkLen = 100;
            P.deflation = 1;
            P.filter_chunkSize = 200000;
            if self.bIsPrefiltered
                return
            end
            if isempty(self.prefilterBufferFilename)
                error('No prefilter bufferfilename was specified when this object was created. Set a bufferfilename first before trying to prefilter!');
            end
%             if exist(self.prefilterBufferFilename, 'file')
%                 warning('PrefilterBufferFile already exists. Delete first before retrying to prefilter!');
%                 return
%             end
            T = self.size(1);
            nC = self.size(2);
            
            self.initializeFilter();
            HD = self.filterObject;
            HD.PersistentMemory = 1;
            dims = [T nC];
            maxDims = [T nC];
            h5type = 'H5T_NATIVE_INT16';
            chunkDims = [P.h5chunkLen nC];
            
            X = mysort.h5.matrixCreate(self.prefilterBufferFilename,...
                                       self.prefilterBufferH5path,...
                                       dims, maxDims, h5type, chunkDims, P.deflation);    
                                   
            mysort.util.filtfilt_chunked(self, X, HD, P.filter_chunkSize, ...
                               self.prefilterScaleFactor);
            % destroy writable object and open in read only mode
            clear X
            readOnly = 1;
            X = mysort.h5.matrix(self.prefilterBufferFilename,...
                                       self.prefilterBufferH5path, readOnly);
            self.prefilterBufferH5matrix = X;
            self.bIsPrefiltered = true;            
        end
        %------------------------------------------------------------------
        function b = isPrefiltered(self)
            b = self.bIsPrefiltered;
        end
        %------------------------------------------------------------------
        function X = getFilteredData(self, varargin)
            if self.bIsPrefiltered
                X = self.getFilterDataFromBufferFile(varargin{:});
            else
                % if not prefiltered, use super method and filter online
                X = getFilteredData@mysort.ds.FilteredDataSourceInterface(self, varargin{:});
            end
        end
        %------------------------------------------------------------------
        function X = getFilterDataFromBufferFile(self, timeIndex, channelIndex)
            X = self.prefilterBufferH5matrix(timeIndex, channelIndex);
            X = double(X)/self.prefilterScaleFactor;
        end        
    end
end