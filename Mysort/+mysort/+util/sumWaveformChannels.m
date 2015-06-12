
function SW = sumWaveformChannels(W, nC)
    warning('This function is depricated! Use mysort.wf.* instead!');
    % Sums the channels of a multichannel waveform
    Tf = size(W,2)/nC;
    SW = zeros(size(W,1), Tf);
    for c=1:nC
        idx = (c-1)*Tf+1:c*Tf;
        SW = SW + W(:,idx);
    end
end