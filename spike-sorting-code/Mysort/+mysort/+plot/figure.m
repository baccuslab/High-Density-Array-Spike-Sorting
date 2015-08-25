
function fh_ = figure(varargin)
    P.color = 'w';
    P.position = [];
    P.width = 600;
    P.height = 300;
    P.w = [];
    P.h = [];
    nargin_ = nargin;
    if nargin_ > 0 && isnumeric(varargin{1}) && length(varargin{1}) == 2
        % special case where it is called as figure([w h], ...)
        P.w = varargin{1}(1);
        P.h = varargin{1}(2);
        varargin = varargin(2:end);
        nargin_ = nargin_ -1;
    end
    
    if nargin_==1
        fh = varargin{1};
        figure(fh);
    else
        [P restargs] = mysort.util.parseInputs(P, varargin, 'split');
        restargs = mysort.util.deflateP(restargs);
        fh = figure(restargs{:});
    end
    if ~isempty(P.h); P.height = P.h; end
    if ~isempty(P.w); P.width  = P.w; end
    set(fh,'color', P.color);
    if isempty(P.position)
        screen_size = get(0, 'ScreenSize');
        position = get(fh, 'position');
        position(1) = round((screen_size(3) - P.width)/2);
        position(3) = P.width;
        position(2) = round((screen_size(4) - P.height)/2);
        position(4) = P.height;        
        set(fh, 'position', position);
    end
    mysort.plot.dataCursorMode(fh);
    if nargout == 1
        fh_ = fh;
    end