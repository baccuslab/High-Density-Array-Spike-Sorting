
function [y BIC minK] = kmeansBIC(X, varargin)
    P.kmin = 1;
    P.kmax = 6;
    P.repeats = 5;
    P.debugFun = [];
    P = mysort.util.parseInputs(P, 'kmeansBIC', varargin);
    
    [N, dim] = size(X);
    nIter = length(P.kmin:P.kmax);

    Y = zeros(N,nIter);
    BIC = zeros(nIter,1);
    idx = 0;
    if ~isempty(P.debugFun); P.debugFun('Clustering k='); end
    opts = statset('Display','off');
        
    for k=P.kmin:min(P.kmax,N)
        if ~isempty(P.debugFun); P.debugFun(sprintf(' %d', k)); end
        idx = idx+1;
        try
            [Y(:,idx), ctrs] = kmeans(X,k,...
                'Distance','sqEuclidean',...
                'Replicates', P.repeats,...
                'Options',opts);
            % Build helper Matrix to construct MU matrix that contains
            % the corresponding cluster center for every point
            M = zeros(N,k);
            for cc = 1:k
                M(Y(:,idx)==cc,cc) = 1;
            end
            MU = M * ctrs;
            NORMS = sqrt(sum(MU.*MU, 2));
            ERR = 1/N * sum(sqrt(sum(((X - MU).*(X - MU)), 2)./(.2*NORMS)));
            BIC(idx) = N*log(ERR) + k*log(N);
        catch
            warning('Empty cluster in all repetitions for k=%d !', k);
            BIC(idx) = inf;
        end
    end
    if ~isempty(P.debugFun); P.debugFun('\n'); end
        
    [ma minK] = min(BIC);

    y = Y(:,minK)+1;
%     minK = P.kmin+floor((minK-1)/P.repeats); 
    
    uy = unique(y);
    for i=1:length(uy)
        y(y==uy(i)) = i;
    end
end


