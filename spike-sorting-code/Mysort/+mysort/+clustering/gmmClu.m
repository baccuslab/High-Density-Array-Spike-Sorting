
function [y BIC minK mu] = gmmClu(X, varargin)
    P.EXPSIG = 6;
    P.kmin = 1;
    P.kmax = 6;
    P.repeats = 3;
    P.debugFun = [];
    P = mysort.util.parseInputs(P, 'gmmClu', varargin);
    
    [N, dim] = size(X);
    nIter = length(P.kmin:P.kmax);

    Y = zeros(N,nIter*P.repeats);
    BIC = zeros(nIter*P.repeats,1);
    idx = 0;
    if ~isempty(P.debugFun); P.debugFun('Clustering k='); end
        
    for k=P.kmin:min(P.kmax,N)
	if ~isempty(P.debugFun); P.debugFun(sprintf(' %d', k)); end
        for r=1:P.repeats
            idx = idx+1;
            %     ini_mu = EXPSIG*randn(k,dim);
            centerIdx = randperm(N);
            ini_mu = X(centerIdx(1:k),:);
            ini_sig = P.EXPSIG*eye(dim);
            [mu{idx}, Sigma, p, LL, Z, Y(:,idx), HCP, HCD] = ...
                mysort.clustering.mixgauss(X,k,'graph3d',0, 'muinit',ini_mu,'siginit',ini_sig);

            %Y(HCP<0.1,idx) = 0;
            %BIC(idx) = -2*L(idx)+k*log(N);
            BIC(idx) = N*log(HCD/N)+5*k*log(N);  % Justify the 5 !!!
        end
    end
    if ~isempty(P.debugFun); P.debugFun('\n'); end
        
    [ma minK] = min(BIC);
%     minK = 9;
    y = Y(:,minK)+1;
    mu = mu{minK};    
    minK = P.kmin+floor((minK-1)/P.repeats); 
end


