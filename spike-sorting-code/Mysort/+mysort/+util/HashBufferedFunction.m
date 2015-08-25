
classdef HashBufferedFunction < mysort.util.CacheableClass
    properties
        hashtable
        fHandle
    end
    methods
        %%% ------------------------------------------------------
        function self = HashBufferedFunction(fHandle, varargin)
            self = self@mysort.util.CacheableClass(varargin{:});
            self.hashtable = java.util.Hashtable();
            self.fHandle = fHandle;
            if self.was_restored; return; end
        end
        
        %%% ------------------------------------------------------
        function init(self)
        end
        
        %%% ------------------------------------------------------
        function R = call(self, hash, varargin)
            if self.hashtable.containsKey(hash)
                R = self.hashtable.get(hash);
                if ischar(R) && strcmp(R, 'empty')
                    R = [];
                end
                if size(R,2) == 1
                    R = R';
                end
                return
            end

            R = self.fHandle(varargin{:});
            if ~isempty(R)
                self.hashtable.put(hash, R);
            else
                self.hashtable.put(hash, 'empty');
            end  
        end
    end
end