
function [v nC2] = v2c(v, nC, c)
    warning('This function is depricated! Use mysort.wf.* instead!');
% select a subset of channels "c" from multichannel waveforms that are
% represented as rows in "v" with "nC" channels concatinated
    nC2 = length(c);
    
    Tf = size(v,2)/nC;
    
    idx = zeros(1,Tf*nC2);
    for i=1:nC2
        idx(1, (i-1)*Tf+1:i*Tf) = (c(i)-1)*Tf+1:c(i)*Tf;
    end
    
    v = v(:, idx);