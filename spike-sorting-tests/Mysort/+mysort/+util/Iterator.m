classdef Iterator < handle
    % This is a standard iterator class
    properties             
        items   % stores the items
        idx     % points to the last returned item
    end
    
    methods
        %------------------------------------------------------------------
        function self = Iterator(items)
            self.items = items;
            self.idx = 0;
        end
        
        %------------------------------------------------------------------
        function b = hasNext(self)
            b = self.idx < length(self.items);
        end        
        
        %------------------------------------------------------------------
        function item = next(self)
            self.idx = self.idx+1;
            if iscell(self.items)
                item = self.items{self.idx};
            else
                item = self.items(self.idx);
            end
        end
        %------------------------------------------------------------------
        function n = getNItems(self)
            n = length(self.items);
        end
        %------------------------------------------------------------------
        function n = length(self)
            n = length(self.items);
        end       
        
    end
end