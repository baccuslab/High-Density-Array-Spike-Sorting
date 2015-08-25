classdef PersistentSpikeSortingContainer < mysort.spiketrain.SpikeSortingContainer
    properties
        fname
        h5path
        gdfIsLoaded
    end
    
    
    methods
        %------------------------------------------------------------------
        function self = PersistentSpikeSortingContainer(fname, h5path, name, varargin)
            if ~exist(fname, 'file')
                [pathstr, fn, ext] = fileparts(fname);
                if ~exist(pathstr, 'dir')
                    mkdir(pathstr);
                end
                %S = mysort.spiketrain.SpikeSortingContainer.emptyStruct();
                S = name;
            else
                S.name = mysort.h5.recursiveLoad(fname, [h5path '/name']);
                S.gdf = [];
                S.details = mysort.h5.recursiveLoad(fname, [h5path '/details']);
            end
            self = self@mysort.spiketrain.SpikeSortingContainer(S, varargin{:});
            self.fname = fname;
            self.h5path = h5path;
            self.gdfIsLoaded = false;
        end
        %------------------------------------------------------------------
        function loadGdf(self)
            self.gdf = double(hdf5read(self.fname, [self.h5path '/gdf'], 'V71Dimensions', true));
            self.gdfIsLoaded = true;
        end
        %------------------------------------------------------------------
        function save2File(self, fname, h5path)
            if nargin == 1
                fname = self.fname;
                h5path = self.h5path;
            end
            save2File@mysort.spiketrain.SpikeSortingContainer(self, fname, h5path);
        end
        
        %------------------------------------------------------------------
        function gdf = getGdf4UnitIdx(self, varargin)
            if ~self.gdfIsLoaded
                self.loadGdf();
            end

            gdf = getGdf4UnitIdx@mysort.spiketrain.SpikeSortingContainer(self, varargin{:});
        end
        %------------------------------------------------------------------
        function T = computeTemplateWfs(self, varargin)
            T = computeTemplateWfs@mysort.spiketrain.SpikeSortingContainer(self, varargin);
            self.save2File();
        end
        %------------------------------------------------------------------
        function T = computeTemplates(self, varargin)
            warning('This function is not properly implemented yet! Arguments are ignored!');
            T = self.computeTemplateWfs();
        end
    end
end