function [P M] = vProjectOnMeanChannelWise(V, nC)
    % Projects the waveforms in rows of V that are concatinated channel-wise
    % onto their mean waveform with respect to channel
    Tf = size(V,2)/nC;
    m = mean(V, 1);
    % build channel wise mean with zeros on the other channels
    M = zeros(nC, Tf*nC);
    for c=1:nC
        idx = (c-1)*Tf+1:c*Tf;
        M(c, idx) = m(idx);
    end
    
    % lazy project. might be worth optimizing a bit
    P = V*M';