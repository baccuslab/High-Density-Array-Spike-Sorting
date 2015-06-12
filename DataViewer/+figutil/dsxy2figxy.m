function varargout = dsxy2figxy(varargin)
    if length(varargin{1})== 1 && ishandle(varargin{1}) && ...
      strcmp(get(varargin{1},'type'),'axes') 
     hAx = varargin{1};
     varargin = varargin(2:end);
    else
     hAx = gca;
    end;
    % Parse either a position vector or two 2-D point tuples
    if length(varargin)==1 % Must be a 4-element POS vector
     pos = varargin{1};
    else
     [x,y] = deal(varargin{:});  % Two tuples (start & end points)
    end
    %% Get limits
    axun = get(hAx,'Units');
    set(hAx,'Units','normalized');  % Need normaized units to do the xform
    axpos = get(hAx,'Position');
    axlim = axis(hAx);              % Get the axis limits [xlim ylim (zlim)]
    axwidth = diff(axlim(1:2));
    axheight = diff(axlim(3:4));
    %% Transform data from figure space to data space
    if exist('x','var')     % Transform a and return pair of points
     varargout{1} = (x-axlim(1))*axpos(3)/axwidth + axpos(1);
     varargout{2} = (y-axlim(3))*axpos(4)/axheight + axpos(2);
    else                    % Transform and return a position rectangle
     pos(1) = (pos(1)-axlim(1))/axwidth*axpos(3) + axpos(1);
     pos(2) = (pos(2)-axlim(3))/axheight*axpos(4) + axpos(2);
     pos(3) = pos(3)*axpos(3)/axwidth;
     pos(4) = pos(4)*axpos(4)/axheight;
     varargout{1} = pos;
    end
    %% Restore axes units
    set(hAx,'Units',axun)