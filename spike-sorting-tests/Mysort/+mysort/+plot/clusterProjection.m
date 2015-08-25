
function [ ] = clusterProjection(spikesX, classes, templates, C, varargin)
    P.name = '';
    P.figure = 1;
    P.figureHandle = [];
    P.axesHandles = [];
    P.plotOnlyDiag = 0;
    P.plotTitleIDs = 0;
    P.plotTitleDiff = 0;
    P.plotTitleRsq = 0;
    P.IDs = [];
    P.iU = [];
    P.colors = {};
    P.binning  = -40:1:40;
    P.xbinning = -40:1:40;
    
    P = mysort.util.parseInputs(P,varargin, 'error');     
    if isempty(P.IDs)
        P.IDs = unique(classes);
    end
    nT = length(P.IDs);
    if nT <= 1
        return
    end
    if nT > 12
        disp('Warning, too many clusters to plot Projections!')
        return
    end
    if nargin <= 3 || isempty(C)
        P.iU = eye(size(spikesX,2));
        pwSpikesX = spikesX;
    else
        if isempty(P.iU)
            U = chol(C);
            P.iU = inv(U);
        end

        pwSpikesX = spikesX*P.iU;
    end
    
    if nargin < 3 || isempty(templates)
        templates = mysort.util.calculateClassMeans(spikesX, classes);
    end
    
    
    if isempty(P.axesHandles)
        if ~isempty(P.figureHandle)
            fh = P.figureHandle;
        elseif P.figure
            fh = mysort.plot.figure('name', 'Cluster Projections');
        end  
        ah = mysort.plot.subplot((nT-1)*(nT-1), 'upperTriangle', 1, 'offsetY', .1);
    else
        ah = P.axesHandles;
    end
    
    if isempty(P.colors)
        P.colors = cell(nT,1);
        for i=1:nT
            if ~isempty(P.IDs)
                P.colors{i} = mysort.plot.vectorColor(P.IDs(i)); 
            else
                P.colors{i} = mysort.plot.vectorColor(i);
            end
        end
    end
    
    normfactor = P.binning(2)-P.binning(1);
    
    i = 1;
    for t1 = 1:nT
        myTID1 = P.IDs(t1);        
        myPwSpikes1 = pwSpikesX(classes == myTID1,:);
        if ~isempty(P.iU)
            myPwTemplate1 = templates(t1,:)*P.iU;
        else
            myPwTemplate1 = templates(t1,:);
        end 
        
        if P.plotOnlyDiag == 0
            t2Set = t1+1:nT;
        else
            t2Set = t1+1;
            t2Set = t2Set(t2Set<=nT);
        end        
        for t2 = t2Set
            myTID2 = P.IDs(t2);        
            myPwSpikes2 = pwSpikesX(classes == myTID2,:);
            if ~isempty(P.iU)
                myPwTemplate2 = templates(t2,:)*P.iU;
            else
                myPwTemplate2 = templates(t2,:);
            end
            Proj = myPwTemplate2-myPwTemplate1;
            projNorm = norm(Proj);
            assert(projNorm > 0, 'Two templates are absolutely identical!')
            Proj = Proj'./projNorm;
            myPSpikes1 = myPwSpikes1*Proj;
            myPSpikes2 = myPwSpikes2*Proj;
                        
            h1 = hist(myPSpikes1, P.binning); h1=h1/sum(h1);
            h2 = hist(myPSpikes2, P.binning); h2=h2/sum(h2);
            bStart = P.binning(find(h1+h2>0,1))-1;
            bEnd   = P.binning(length(P.binning) - find(fliplr(h1+h2)>0,1) + 1)+1;
            
            m1 = mean(myPSpikes1);
            m2 = mean(myPSpikes2);
            std1 = 1;%std(myPSpikes1);
            std2 = 1;%std(myPSpikes2);
            D = abs(m1-m2);
            
            %Rsquare = 1 - (SSerr/SStot);
            f1 = normpdf_inline(P.binning,m1,std1);
            ssErr1 = sum(( h1-f1       ).^2);
            ssTot1 = sum(( h1-mean(h1) ).^2);
            rsq1 = 1 - ssErr1/ssTot1;
            
            f2 = normpdf_inline(P.binning,m2,std2);
            ssErr2 = sum(( h2-f2       ).^2);
            ssTot2 = sum(( h2-mean(h2) ).^2);
            rsq2 = 1 - ssErr2/ssTot2;
            
            set(ah(i), 'nextplot', 'add');
            bar(ah(i), P.binning, h1); bar(ah(i),P.binning, h2);
            axis(ah(i), [bStart bEnd 0 .45*normfactor]);
            
            h = findobj(ah(i),'Type','patch');
            % Warning, the ordering of bar plots and pathes in "h" is
            % reversed!! 
            set(h(2),'FaceColor', P.colors{t1})%,'EdgeColor','w')   
            set(h(1),'FaceColor', P.colors{t2})%,'EdgeColor','w')               
            
            plot(ah(i), P.xbinning, normfactor*normpdf_inline(P.xbinning,m1,std1),'-k','lineWidth',1.5);
            plot(ah(i), P.xbinning, normfactor*normpdf_inline(P.xbinning,m2,std2),'-k','lineWidth',1.5);     
            titleStr = '';
            set(ah(i), 'yticklabel', []);
            if P.plotTitleIDs
                titleStr = sprintf('%d|%d', myTID1, myTID2);
            end
            if P.plotTitleDiff
                titleStr = [titleStr sprintf('|%.1f', D)];
            end
            if P.plotTitleRsq
                titleStr = [titleStr sprintf(' (%.1f|%.1f)', rsq1, rsq2)];
            end
            if ~isempty(titleStr)
                title(ah(i), titleStr);
            end
            i = i+1;
        end
    end
    
    %----------------------------------------------------------------------
    function y = normpdf_inline(x, mu, sigma)
        y = exp(-.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
    end
end