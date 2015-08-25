function P = mds2D(X, T, varargin)
    P.nNeighbors = 3;
    P.fh = [];
    P.ah = [];
    P = mysort.util.parseInputs(P, varargin);
    
    if isempty(P.ah) && isempty(P.fh)
        P.fh = mysort.plot.figure('w', 1100, 'h', 800);
    end

    if isempty(P.ah) 
        P.ah = axes();
    end
    
    K = P.nNeighbors;
    
    DT = pdist(T);
    YT = mdscale(DT,2, 'criterion', 'metricstress');
    
    nT = size(T,1);
    
    N = zeros(nT, K);
    for i=1:nT
        di = sum( (repmat(T(i,:), nT,1) - T).^2 , 2);
        [sdi, sidx] = sort(di, 'ascend');
        N(i,:) = sidx(2:K+1);
    end
    
    
    plot(P.ah, YT(:,1), YT(:,2), 'kx', 'linewidth',3)
    set(P.ah, 'nextplot','add');
    
    for i=1:nT
        for k=1:K
            plot([YT(i,1) YT(N(i,k),1)], [YT(i,2) YT(N(i,k),2)], 'r-')
        end
    end

    
    
    
    
    