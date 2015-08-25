
function aH = XIvsFAxesHandles(nXI, nF, varargin)
    P.offsetX = .015;
    P.offsetY = .08;
    P.spacerX = .025;
    P.spacerY = .1;
    P.bigSpacerX = .05;
    P.bigSpacerY = .15; 
    P = mysort.util.parseInputs(P,'plotXIvsFAxesHandles',varargin);
    
    w = (1-P.offsetX-P.bigSpacerX-(nXI)*P.spacerX)/(nXI+1);
    h = max((1-P.offsetY-P.bigSpacerY-(nF)*P.spacerY)/(nF+1), .01); 

    % calculate XCoords and YCoords:
    for i=1:nXI+1
        if i==1
            xC(i) = P.offsetX;
        else
            xC(i) = P.offsetX + P.bigSpacerX + (i-2)*P.spacerX+(i-1)*w;
        end
    end
    for i=1:nF+1
        if i==1
            yC(i) = P.offsetY + P.bigSpacerY + (nF-1)*P.spacerY+nF*h;
        else
            yC(i) = yC(1)-P.bigSpacerY - (i-1)*(h)-(i-2)*P.spacerY;
        end
    end
    aH = [];
    count = 1;
    for j=1:nF+1
        for i=1:nXI+1
            if i+j>2
                aH(count) = axes('Units','normalized','position',[xC(i) yC(j) w h]);
                count = count +1;
            end
        end
    end
end