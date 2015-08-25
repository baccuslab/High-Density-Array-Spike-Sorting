function vidx = vSubIdx(V_or_Tf, nC, idx)
    % Computes the subidx vidx for all channels in V (concatenated) when
    % idx is just the index into one channel.
    if size(V_or_Tf,1) == 1 && size(V_or_Tf,2) == 1
        Tf = V_or_Tf;
    else
        Tf = size(V_or_Tf,2)/nC;
        assert(round(Tf)==Tf, 'dim2 of V does not match number of channels!');
    end
    nI = length(idx);
    vidx = zeros(1, nI*nC);
    
    for c=1:nC
        vidx((c-1)*nI+1:nI*c) = (c-1)*Tf + idx;
    end