classdef FileBasedMatrix < mysort.ds.Matrix
                           
                    
    properties
        fname
        fpath
    end
    methods
        %------------------------------------------------------------------
        function self = FileBasedMatrix(ffname, X, samplesPerSecond, varargin)
%             P.hpf = 300;
%             P.lpf = 6000;
%             P.filterOrder = 6;
%             P.filterType = 'butter';
%             [path, fname, extension] = fileparts(ffname);
%             prefilterBufferFilename = fullfile(path, fname, 'prefilter', ...
%                 [num2str(P.filterOrder) '_' num2str(P.hpf)...
%                 '_' num2str(P.lpf) '.h5']);
%             prefilterBufferH5path = '/prefilterdData';
%             
%             if ischar(samplesPerSecond)
%                 samplesPerSecond = hdf5read(ffname, samplesPerSecond);
%             end
%             
%             filterFactory = mysort.mea.FilterWrapper(P.hpf, ...
%                 P.lpf, samplesPerSecond, P.filterOrder, P.filterType);
            self = self@mysort.ds.Matrix(X, samplesPerSecond, varargin{:});
            self.fname = ffname;

            assert(ndims(X)==2, 'data must be either a ds.Interface or a matrix!');            
        end
%         %------------------------------------------------------------------
%         function X = getUnfilteredData(self, varargin)
%             X = self.X.getData(varargin{:});
%         end
    end
end