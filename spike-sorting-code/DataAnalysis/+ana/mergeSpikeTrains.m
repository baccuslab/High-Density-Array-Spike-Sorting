function [groups maxT D] = mergeSpikeTrains(gdf, T, varargin)
    P.binSize = 4;
    P.maxLag = 4;
    P.nMaxSpikesPerSpikeTrain = 2000;
    P.minPercentOfEqualSPikes = 20;
    P = mysort.util.parseInputs(P, varargin, 'error');
    absMaxPeaksOfRemainingTemplates = max(abs(mysort.wf.t2v(T)), [], 2);
    
    classes = unique(gdf(:,1));
    nT = length(classes);
    XC = zeros(nT,nT);
    keepTemplates = ones(1,nT);
    ci = 1:nT;
    for i=1:nT
        if ~keepTemplates(i)
            continue
        end
        s1 = gdf(gdf(:,1)==classes(i),2);
        if length(s1)>P.nMaxSpikesPerSpikeTrain
            s1 = s1(1:P.nMaxSpikesPerSpikeTrain);
        end
        for j=i+1:nT
            if ~keepTemplates(j)
                continue
            end
            s2 = gdf(gdf(:,1)==classes(j),2);
            if length(s2)>P.nMaxSpikesPerSpikeTrain
                s2 = s2(1:P.nMaxSpikesPerSpikeTrain);
            end
            [xc, E_xc, var_xc, conf_95_xc, bin_centers, edges, xc_unnormalized] = ...
                mysort.spiketrain.xcorr(s1, s2, 'maxLag', P.maxLag, 'binSize', P.binSize);
            XC(i,j) = max(100*xc_unnormalized/length(s1));
            XC(j,i) = max(100*xc_unnormalized/length(s2));
            
            if XC(i,j) > P.minPercentOfEqualSPikes ||...
               XC(j,i) > P.minPercentOfEqualSPikes
                if absMaxPeaksOfRemainingTemplates(i) > absMaxPeaksOfRemainingTemplates(j)
                    keepTemplates(j) = 0;
                    ci(ci==j) = i;
                else
                    keepTemplates(i) = 0;
                    ci(ci==i) = j;
                    break
                end
                % mysort.plot.templates2D(T(:,:,[i j]), MES.electrodePositions, 15, [], 'stacked', 0);
            end
        end
    end   
    uGroups = unique(ci);
    groups = cell(1, length(uGroups));
    maxT = zeros(1,length(uGroups));
    for i=1:length(uGroups)
        groups{i} = find(ci==uGroups(i));
        maxT(i) = uGroups(i);
    end
    
    if nargout > 2
        D.XC = XC;
        D.keepTemplates = keepTemplates;
        D.ci = ci;
    end