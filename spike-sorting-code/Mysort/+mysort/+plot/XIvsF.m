
function [RES] = XIvsF(XI,F, varargin)
    % XIvsF: 1st Dim = samples, 2nd Dim = f, 3rd Dim = xi
    % the f-dimension are the new "channels". The old waveform xi is now
    % represented in the filter space as an #f channel signal
    P.figure = 1;
    P.figureHandle= -1;
    P.XIvsF = [];
    P.color = '-b';
    P.diagonalColor = [];
    P.holdOn = 0;
    P.mcSpacerXI = -1;
    P.mcSpacerF = -1;  
    P.ymax = 1.1;
    P.ymin = -0.5;
    P.axesHandles = [];
    P.sampleInterval = 1;
    P.plotXIF=1;   
    P.axistight = 1;
    P.title = 0;
    P.nC = 1;
    import mysort.*
    P = util.parseInputs(P,'plotXIvsF',varargin);
    
    RES.fh = 0;
    if P.figureHandle > -1
        RES.fh = P.figureHandle;
    elseif P.figure
        RES.fh = mysort.plot.figure('color','w');
    end
    
    if ndims(XI) == 2
        XI = util.m2t(XI, P.nC);
    end
    
    if ndims(F) == 2
        F = util.m2t(F, P.nC);
    end
    
    
    if isempty(P.XIvsF)
        P.XIvsF = util.calculateXIvsF(util.t2m(XI),util.t2m(F), P.nC, 1); 
    end

    Tf = size(F, 1);
    nC = size(F, 2);
    nF = size(F, 3);
    nT = size(XI, 3);
    maxXI = max(XI(:));
    minXI = min(XI(:));
    maxF = max(F(:));
    minF = min(F(:));
    
    if isempty(P.axesHandles)
        if P.title==1
            P.axesHandles = mysort.plot.XIvsFAxesHandles(nT, nF,'plotF',P.plotXIF);
        else
            P.axesHandles = mysort.plot.XIvsFAxesHandles(nT, nF,'plotF',P.plotXIF,...
                'spacerX',.005, 'spacerY', .005, ...
                'bigSpacerY', .05);            
        end
    end

    if nargout>0
        RES.XIvsF = P.XIvsF;
        RES.axesHandles = P.axesHandles;
    end
    
    xdata = (1:Tf)*P.sampleInterval;
    
    for x=1:nT
        a=P.axesHandles( x); axes(a);
        if P.holdOn ==1; hold on; end
        P.mcSpacerXI = plot.mc(XI(:,:,x)','figure',0,'color',{P.color},'spacer',P.mcSpacerXI,'showZeroLine',0,'srate',P.sampleInterval);
        if P.title
            title( ['$$\xi^{' num2str(x) '}$$'],'Interpreter','latex','FontSize',14)
        end
%         if f<nF
%             set(a,'xticklabel', []);
%         end
        if x>1
            set(a,'yticklabel',[]);
        end
        axis([0 P.sampleInterval*Tf minXI nC*P.mcSpacerXI+maxXI]);
             
        for f=1:nF
            filterPlotIdx = nT+1 + (f-1)*(nT+1);
            if x == 1
                % Plot Filter on the left side
                a=P.axesHandles(filterPlotIdx); axes(a);
                if P.holdOn ==1; hold on; end
                P.mcSpacerF = plot.mc(F(:,:,f)','figure',0,'color',{P.color},'spacer',P.mcSpacerF,'showZeroLine',0,'srate',P.sampleInterval);
                if P.title
                    title( ['$$f^{' num2str(f) '}$$'],'Interpreter','latex','FontSize',14)
                end
                if f~=nF
                    set(a,'xticklabel', []);
                end
                if f>1
                    %set(a,'yticklabel', []);
                end
                if nC>1; ymax = nC*P.mcSpacerF+maxF; else ymax = maxF; end
                axis([0 P.sampleInterval*Tf minF ymax]);     
            end
            
            pIdx = filterPlotIdx + x;
            a=P.axesHandles(pIdx); axes(a);
            if P.holdOn ==1; hold on; end
            tmp = squeeze(P.XIvsF(:,f,x));

            if ~isempty(P.diagonalColor) && (mod(f-1,nT)==x-1)
                pColor = P.diagonalColor;
            else
                pColor = P.color;
            end
            plot([-(Tf-1):(Tf-1)], tmp, pColor);
            if P.title
                title( ['$$(\xi^{' num2str(x) '}\star f^{' num2str(f) '})_{\tau}$$'],'Interpreter','latex','FontSize',14)
            end
%             if f==1
%                 title(sprintf('F: %d', x));
%             end
            if f~=nF
                set(a,'xticklabel', []);
            end
            if P.axistight
                axis tight
            else
                if x==1
    %                 ylabel(sprintf('XI: %d', f));
                else
                    set(a,'yticklabel', []);
                end
                axis([-(Tf-1) (Tf-1) P.ymin P.ymax]);
            end
        end
    end
end
