function [mi, ma, mi_idx, ma_idx] = tMinMaxPerTemplate(T)
    % Returns the minima and maxima for each template in T on each channel
    % Input: 
    %    T  - tensor (time x channels x templates)
    % Output:
    %    mi - minima for each template on each channel
    %    ma - maxima
    %    mi_idx - temporal index for that minimum on its channel
    %    ma_idx - same for max
    
    [nS, nC, nT] = size(T);
    mi = zeros(nT, nC);
    ma = zeros(nT, nC);
    mi_idx = zeros(nT,nC);
    ma_idx = zeros(nT,nC);
    
    % use for loops to catch special cases (1 channel, 1 template)
    for c=1:nC
        for t=1:nT
            [mi(t,c), mi_idx(t,c)] = min(T(:,c,t));
            [ma(t,c), ma_idx(t,c)] = max(T(:,c,t));
        end
    end
    