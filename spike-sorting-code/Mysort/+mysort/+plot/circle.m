
function circle(ah, c, r, varargin)
    if ~ishandle(ah)
        varargin = [r varargin];
        r = c;
        c = ah;
        ah = gca;
    end
    if nargin < 1
        c = [0 0];
    end
    if nargin < 2
        r = 1;
    end
    th = 0:pi/50:2*pi;
    xunit = r*cos(th)+c(1);
    yunit = r*sin(th)+c(2);
    plot(ah, xunit, yunit, varargin{:});
end