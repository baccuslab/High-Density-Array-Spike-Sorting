function cidx = vSubChannelIdx(V_or_Tf, nC, channel_idx)
    % Computes the subidx cidx for the complete timerange for a subset of
    % all channels in V (concatenated) when
    % channel_idx is just the channel index
    if size(V_or_Tf,1) == 1 && size(V_or_Tf,2) == 1
        Tf = V_or_Tf;
    else
        Tf = size(V_or_Tf,2)/nC;
        assert(round(Tf)==Tf, 'dim2 of V does not match number of channels!');
    end
    if nargin == 2
        channel_idx = 1:nC;
    end
    nCout = length(channel_idx);
    cidx = zeros(1, Tf*nCout);
    timeIndex = 1:Tf;
    for c=1:nCout
        c1 = (c-1)*Tf+1;
        c2 = c1+Tf-1;
        cidx(c1:c2) = (channel_idx(c)-1)*Tf + timeIndex;
    end