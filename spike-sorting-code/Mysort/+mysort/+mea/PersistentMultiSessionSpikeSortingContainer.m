classdef PersistentMultiSessionSpikeSortingContainer < mysort.mea.MultiSessionSpikeSortingContainer
    properties
        filename
    end
    
    methods
        %------------------------------------------------------------------
        function self = PersistentMultiSessionSpikeSortingContainer(filename, varargin)
            if nargin == 1
                varargin = {filename};
            else
                [pathstr, name, ext] = fileparts(filename);
                if ~exist(pathstr, 'file')
                    try
                        mkdir(pathstr);
                    end            
                end
                if exist(filename, 'file')
                    S = mysort.mea.PersistentMultiSessionSpikeSortingContainer.loadSFromFile_(filename);
                    varargin = {S};
                end
            end
            self = self@mysort.mea.MultiSessionSpikeSortingContainer(varargin{:});
            self.filename = filename;
        end
        %------------------------------------------------------------------
        function laodFromFile(self)
%             S = mysort.h5.recursiveLoad(self.filename, self.h5path);
            S = self.loadSFromFile_(self.filename);
            self.fromStruct(S);
        end
       
        %------------------------------------------------------------------
        function save2File(self, filename)
            S = self.toStruct();
            %mysort.h5.recursiveSave(filename, S, h5path);
            save(filename, 'S');
        end
        %------------------------------------------------------------------
        function save(self)
            self.save2File(self.filename);
        end
    end
    methods (Static) 
        %------------------------------------------------------------------
        function S = loadSFromFile_(filename)
%             S = mysort.h5.recursiveLoad(self.filename, self.h5path);
            S = load(filename, '-mat', 'S');
            S = S.S;
        end 
    end
end