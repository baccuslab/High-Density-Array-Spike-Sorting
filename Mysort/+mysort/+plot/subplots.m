
function ax = subplots(nP, varargin)
    P.offsetX = .06;
    P.offsetY = .06;
    P.spacerTop   = .08;
    P.spacerRight = .0;
    P.spacerX = .05;
    P.spacerY = .08;
    P.labels  = {};
    P.preferVertical = 1;
    P.figureTitle = [];
    P.matrix = 0;
    P.holdOn = 1;
    P.upperTriangle = false;
    P = mysort.util.parseInputs(P,varargin,'error');
    
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
    assert(nR < 37 && nC <37, 'Too many rows or columns to plot!');
    
    w = max(0.01, (.99-P.offsetX-(nC-1)*P.spacerX-P.spacerRight)/(nC));
    h = max(0.01, (.99-P.offsetY-(nR-1)*P.spacerY-P.spacerTop)/(nR));
    
    % check if spacers are too big
    maxHeight = P.offsetY + (nR-1)*(P.spacerY + h);
    if maxHeight > 1
        P.spacerY = 0;
        h = max(0.01, (1-2*P.offsetY)/(nR));
    end
    maxWidth = P.offsetX + (nC-1)*(P.spacerX + w);
    if maxWidth > 1 || P.spacerX > w
        P.spacerX = 0;
        w = max(0.01, (1-1.5*P.offsetX-(nC-1)*P.spacerX)/(nC));
    end
    
	if ~isempty(P.figureTitle)
        % Correct hight, to leave room for title
        h = max(0.01, h - .05/nR);
        mysort.plot.figureTitle(P.figureTitle);
	end
    % calculate XCoords and YCoords:
    count = 1; xC = zeros(nC,1); yC = zeros(nR,1); ax = zeros(nP,1);
    for j=nR:-1:1
        for i=1:nC
            if count > nP
                continue
            end            
            if ~P.upperTriangle || i>=(nR-j+1)
                xC(i) = P.offsetX + (i-1)*(P.spacerX + w);
                yC(j) = P.offsetY + (j-1)*(P.spacerY + h);
                ax(count) = axes('Units','normalized','position',[xC(i) yC(j) w h]);
                if P.holdOn
                    set(ax(count), 'nextplot', 'add');
                end
                if ~isempty(P.labels)
                    set(ax(count),'Units', 'pixel');
                    pos = get(ax(count), 'Position' );
                    set(ax(count),'Units', 'normalized');
                    ann = annotation('textbox',[0 .1 0 .1],...
                        'Units','pixel',...
                        'FitBoxToText','off',...
                        'String',{['\textbf{' P.labels{count} '}']},...
                        'FontSize', 14,...
                        'Interpreter', 'latex',... 
                        'LineStyle','none');
                    set(ann, 'Position', [max(0,pos(1)-48) pos(2)+pos(4)+5 30 30]);
                    set(ann, 'Unit', 'normalized');
                end
                count = count +1;
            end
        end
    end
    if P.upperTriangle
        ax = ax(ax>0);
    end
    if P.matrix
        if P.upperTriangle
            tmp = ax;
            ax = nan([nC nR]);
            count = 0;
            for i=1:nC
                for j=i:nR
                    count = count+1;
                    ax(i,j) = tmp(count);
                end
            end
        else
            ax = reshape(ax, [nC nR])';
        end
    end
