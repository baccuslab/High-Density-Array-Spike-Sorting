
function h = verticalLines(ah,xin,mnmx,varargin)

%plot_lines(xin)   Plots vertical lines at points defined on the x-axis
%   using red solid lines and the current y-axis limits using a single line
%   handle
%
%plot_lines(xin,[a b])  Plots the lines from a to b
%
%plot_lines(xin,[],lntype)  Plots the lines using the current y-axis limits
%using the line type deined by lntype
%
%h = plot_lines(xin,[a b],lntype)  Returns the handle to the vertcal lines
%
%
%Ex.)
%figure(1);clf;
%plot(randn(1,1000));
%x = [10 170 300 450 600 800 990];
%hold on;
%plot_lines(x,[],'b:');

%
%Written by S. Wegerich, SmartSignal, Corp. 09/03
%
% Modified by Felix Franke

if ~ishandle(ah)
    if nargin > 2
        varargin = [{mnmx} varargin];
    end
    if nargin > 1
        mnmx = xin;
    end
    xin = ah;
    ah = gca;
end

if((nargin<3)|isempty(mnmx))
    mnmx = get(ah,'YLim');
end

% if((nargin<3)|isempty(lntype))
%     lntype = '-r';
% end

xin = xin(:)';
L = length(xin);
x = reshape([xin;xin;ones(1,L)*nan;],L*3,1);

x = x(1:end-1);
x = x(:);
mnmx = mnmx(:)';

y = repmat([mnmx nan],1,L);
y = y(1:end-1);

h = plot(ah, x,y,varargin{:});