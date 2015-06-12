classdef GuiElement < handle
    properties
        % gui elements
        Parent 
        
        % internals
        bIsVisible
    end
    methods (Abstract)
        makeLayout(self)
        update(self)
    end
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = GuiElement(varargin)
            P = struct();
            P = guiutil.GuiElement.checkInputs(P, varargin);
            assert(length(P.Parent) == 1, 'Only one Parent allowed!');
            bIsAllowedParent = ishandle(P.Parent) || isa(P.Parent, 'uiextras.Container');
            assert(bIsAllowedParent, 'Parent is not an allowed gui handle!');
            self.Parent = P.Parent;
            
            assert(islogical(P.visible), 'visible must be logical!');
            self.bIsVisible = P.visible;
        end
        %------------------------------------------------------------------
        function setVisibility(self, b)
            if b && ~self.bIsVisible
                self.bIsVisible = b;
                self.update();
            else
                self.bIsVisible = b;            
            end
        end      
    end
    %----------------------------------------------------------------------
    methods (Static)
        %------------------------------------------------------------------
        function P = checkInputs(P, V)
            P.Parent = [];
            P.visible = true;            
            if ~isempty(V)
                while ~isempty(V) && ~ischar(V{1})
                    if ~isempty(V{1})
                        if (~isempty(V{1}) && ishandle(V{1})) || ...
                                  isa(V{1}, 'uiextras.Container')
                            P.Parent = V{1};
                            V(1) = [];                          
                        elseif isa(V{1}, 'mysort.ds.DataSourceInterface')
                            P.dataSource = V{1};
                            V(1) = [];
                        else
                            break
                        end
                    end
                end
            else
                P.Parent = mysort.plot.figure();
            end
            P = mysort.util.parseInputs(P, V, 'error');          
        end        
    end
end

