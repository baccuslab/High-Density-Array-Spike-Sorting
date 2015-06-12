function [m idx] = mFindMinInterp(M)
    % find the minumum of the multichannel waveform in M. individual rows
    % of M correspond to individual channels
    m = zeros(size(M,1),1);
    idx = m;
    for i=1:size(M,1)
        [m(i) idx(i)] = mysort.util.findMinInterp(M(i,:));
    end