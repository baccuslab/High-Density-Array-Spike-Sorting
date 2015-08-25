function ellipse(c, r, scx, scy, varargin)
    if nargin < 1
        c = [0 0];
    end
    if nargin < 2
        r = 1;
    end
    th = 0:pi/50:2*pi;
    xunit = r*cos(th)*scx+c(1);
    yunit = r*sin(th)*scy+c(2);
    plot(xunit, yunit, varargin{:});
end