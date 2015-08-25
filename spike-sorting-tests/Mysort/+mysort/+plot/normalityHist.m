function edges = normalityHist(x, varargin)
    P.axesHandles = [];
    P.plotNormFit = 1;
    P.edges = [];
    P.nBins = 20;
    P.oversampleFit = 4;
    P.normFitPara = {'color', 'r', 'linewidth', 2};
    [P uP] = mysort.util.parseInputs(P, varargin, 'split');
    barparams = mysort.util.deflateP(uP);
    if isempty(x)
        return
    end
    assert(~any(isnan(x)), 'nan in x!');
    mi = min(x);
    ma = max(x);
    m = mean(x);
    st = 1; % dont use std(x) !! the idea is to check against normality
    if mi == ma 
        mi = mi-1;
        ma = ma+1;
    end
    if isempty(P.edges)
        P.edges = linspace(mi, ma, P.nBins);
    end
        
    dX = P.edges(2)-P.edges(1);
    edges = [P.edges P.edges(end)+dX];
    middles = (edges(1:end-1)+edges(2:end))/2;
    
    
    n = histc(x,edges);
    N = sum(n(1:end-1));
    if isempty(P.axesHandles)
        P.axesHandles = gca;
    end
    
    bar(P.axesHandles, middles, n(1:end-1), barparams{:});
    if P.plotNormFit
        dXf = dX/P.oversampleFit;
        f_middles = middles(1):dXf:middles(end);
        f = N*dX*normpdf_inline(f_middles,m,st);
        set(P.axesHandles, 'nextplot', 'add');
        plot(P.axesHandles, f_middles, f, P.normFitPara{:})
    end
    axis(P.axesHandles, 'tight');
    set(P.axesHandles, 'xlim', [edges(1) edges(end)])
    
    if nargout == 1
        edges = P.edges;
    end
        
%     ssErr1 = sum(( h1-f1       ).^2);
%     ssTot1 = sum(( h1-mean(h1) ).^2);
%     rsq1 = 1 - ssErr1/ssTot1;
            
            
    %----------------------------------------------------------------------
    function y = normpdf_inline(x, mu, sigma)
        y = exp(-.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
    end
end
