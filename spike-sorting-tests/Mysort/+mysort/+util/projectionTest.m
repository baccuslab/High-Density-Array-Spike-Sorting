function [M STD D merged] = projectionTest(spikesX, gdf, templates, C, varargin)
    P.iU = [];
%     P.binning  = -40:1:40;
%     P.xbinning = -40:1:40;
    
    P = mysort.util.parseInputs(P,varargin, 'error');     
    classes = unique(gdf(:,1));
    nT = length(classes);
    if nT <= 1
        return
    end
    if nargin <= 3 || isempty(C)
        pwSpikesX = spikesX;
        P.iU = eye(size(spikesX,2));
    else
        if isempty(P.iU)
            U = chol(C);
            P.iU = inv(U);
        end
        pwSpikesX = spikesX*P.iU;
    end
    
    if nargin < 3 || isempty(templates)
        templates = mysort.util.calculateClassMeans(spikesX, gdf(:,1));
    end
    
%     normfactor = P.binning(2)-P.binning(1);
    
    D = zeros(nT, nT);
    M = zeros(nT, nT);
    STD = zeros(nT, nT, 2);
    i = 1; 
    for t1 = 1:nT
        myTID1 = classes(t1);        
        myPwSpikes1 = pwSpikesX(gdf(:,1) == myTID1,:);
        if ~isempty(P.iU)
            myPwTemplate1 = templates(t1,:)*P.iU;
        else
            myPwTemplate1 = templates(t1,:);
        end 
        
        t2Set = t1+1:nT;
        for t2 = t2Set
            myTID2 = classes(t2);        
            myPwSpikes2 = pwSpikesX(gdf(:,1) == myTID2,:);
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
                        
%             h1 = hist(myPSpikes1, P.binning); h1=h1/sum(h1);
%             h2 = hist(myPSpikes2, P.binning); h2=h2/sum(h2);
%             bStart = P.binning(find(h1+h2>0,1))-1;
%             bEnd   = P.binning(length(P.binning) - find(fliplr(h1+h2)>0,1) + 1)+1;
            
            m1 = mean(myPSpikes1);
            m2 = mean(myPSpikes2);
            M(t1,t2,1) = m1;
            M(t1,t2,2) = m2;
            STD(t1,t2,1) = std(myPSpikes1);
            STD(t1,t2,2) = std(myPSpikes2);
            D(t1,t2) = abs(m1-m2);
            
%             %Rsquare = 1 - (SSerr/SStot);
%             f1 = normpdf_inline(P.binning,m1,std1);
%             ssErr1 = sum(( h1-f1       ).^2);
%             ssTot1 = sum(( h1-mean(h1) ).^2);
%             rsq1 = 1 - ssErr1/ssTot1;
%             
%             f2 = normpdf_inline(P.binning,m2,std2);
%             ssErr2 = sum(( h2-f2       ).^2);
%             ssTot2 = sum(( h2-mean(h2) ).^2);
%             rsq2 = 1 - ssErr2/ssTot2;
%             
         
            i = i+1;
        end
    end
    

    
%     %----------------------------------------------------------------------
%     function y = normpdf_inline(x, mu, sigma)
%         y = exp(-.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
%     end
end