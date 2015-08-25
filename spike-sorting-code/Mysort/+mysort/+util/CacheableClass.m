
classdef CacheableClass < mysort.util.DebuggableClass
    properties (Constant)

    end
    properties

    end
    
    methods (Abstract)
        init(self)
    end
    
    methods
        %%% ------------------------------------------------------
        function self = CacheableClass(varargin)
            self = self@mysort.util.DebuggableClass(varargin{:});
            if self.was_restored; return; end
            self.P.cache_name = [];
            self.P.cache_path = [];
            self.P.cacheable = false;
            self.P.cache_reset = false;
            suppress_warning = 1;
            self.P = mysort.util.parseInputs(self.P, 'CacheableClass', varargin, suppress_warning);         
            
            if self.P.cache_reset
                try 
                    self.debugout('Resetting cache.', self.LEVEL_FLOW);
                    if exist(get_cache_fullfile(self), 'file')
                        delete();
                    end
                catch
                end
            elseif self.P.cacheable
                cached_self = self.load_from_cache();
                if ~isempty(cached_self)
                    self = cached_self;
                    self.P.cacheable = true;
                    return
                end
            end
            self.init(varargin{:});
        end
        
        %%% ------------------------------------------------------
        function cache_me(self)
            self.cache_check();
            save(self.get_cache_fullfile(), 'self');          
        end

        %%% ------------------------------------------------------
        function me = load_from_cache(self)
            self.cache_check();
            try 
                self.debugout('Trying to load cached class... ', self.LEVEL_FLOW);
                me = load(self.get_cache_fullfile());
                me = me.self;
                self.debugout('done.', self.LEVEL_FLOW);
            catch
                self.debugout('Failed to load cached class.', self.LEVEL_FLOW);
                me = [];
            end                      
        end        
        
        %%% ------------------------------------------------------
        %%% -----------------HELPER FUNCTIONS---------------------
        %%% ------------------------------------------------------
        function f = get_cache_fullfile(self)
            f = fullfile(self.P.cache_path, ['cached_' self.P.cache_name '.mat']);
        end
        
        %%% ------------------------------------------------------
        function cache_check(self)
            assert(self.P.cacheable, 'This class was not set to be cacheable!');
            assert(~isempty(self.P.cache_path), 'To make a class cacheable, provide cache_path!');
            assert(~isempty(self.P.cache_name), 'To make a class cacheable, provide caceh_name!');
        end
    end
end