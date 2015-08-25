
function ax = subplot2(nP, varargin)
    % Subplot control in pixels
    
    % bottom, top, left, right
    P.marginBottom = 70;
    P.marginTop    = 50;
    P.marginLeft   = 70;
    P.marginRight  = 0;
    P.spacerX = 80;
    P.spacerY = 100;
    P.preferVertical = 1;
    P.figureTitle = [];
    P.upperTriangle = false;
    P.labels = {};
    P.fig = [];
    P.max = false;
    P = mysort.util.parseInputs(P,'subplot2',varargin);
    
    if isempty(P.fig)
        P.fig = gcf;
    end
    
    if length(nP) == 2
        nR = nP(1);
        nC = nP(2);
        nP = nR*nC;
    elseif P.preferVertical
        nC = floor(sqrt(nP));
        nR = ceil(nP/nC);
    else
        nR = floor(sqrt(nP));
        nC = ceil(nP/nR);
    end
    
    assert(isempty(P.labels) || length(P.labels)>=nP, 'If labels is provided it must contain one label per subplot!');
    
    if ~isempty(P.figureTitle)
        % Correct hight, to leave room for title
        P.marginTop = P.marginTop + 30;
        mysort.plot.figureTitle(P.fig, P.figureTitle, true);
    end
    if P.max
        P.marginTop = 0;
        P.marginBottom = 0;
        P.marginLeft = 0;
        P.marginRight = 0;
        P.spacerX = 0;
        P.spacerY = 0;
    end
    set(P.fig, 'Units', 'pixel');
    pos = get(P.fig, 'Position');
    W = pos(3);
    H = pos(4);
    w = (W - P.marginLeft - P.marginRight  -(nC-1)*P.spacerX)/(nC);
    h = (H - P.marginTop  - P.marginBottom -(nR-1)*P.spacerY)/(nR);
    w = max(10, w);
    h = max(10, h);
    
    % calculate XCoords and YCoords:
    count = 1; xC = zeros(nC,1); yC = zeros(nR,1); ax = zeros(nP,1);
    for j=nR:-1:1
        for i=1:nC
            if ~P.upperTriangle || i>=(nR-j+1)
                xC(i) = P.marginLeft + (i-1)*(P.spacerX + w);
                yC(j) = P.marginBottom + (j-1)*(P.spacerY + h);
                ax(count) = axes('Units','pixel','position',[xC(i) yC(j) w h]);
                if ~isempty(P.labels)
                    set(ax(count),'Units', 'pixel');
                    pos = get(ax(count), 'Position' );
%                     set(ax(count),'Units', 'normalized');
%                     ann = annotation('textbox',[0 .1 0 .1],...
%                         'Units','pixel',...
%                         'FitBoxToText','off',...
%                         'String',{['\textbf{' P.labels{count} '}']},...
%                         'FontSize', 14,...
%                         'Interpreter', 'latex',... 
%                         'LineStyle','none');
%                     set(ann, 'Position', [max(0,pos(1)-48) pos(2)+pos(4)+5 30 30]);
%                     set(ann, 'Unit', 'normalized');
                end
                count = count +1;
            end
            if count > nP
                break
            end
        end
    end
    if P.upperTriangle
        ax = ax(ax>0);
    end
    mysort.plot.figureChildrenSet(P.fig, 'Units', 'normalized');
    axes(ax(1));