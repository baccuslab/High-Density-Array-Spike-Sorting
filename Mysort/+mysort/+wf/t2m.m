function M = t2m(T)
    % converts the multi channel waveform tensor T into the matrix
    % representation for waveforms where every row of M represents a single
    % channel waveform and multiple channels of one waveform are below each
    % other in M.
    % If nS is the number of waveforms, Tf is the number of samples per
    % individual channel and nC is the number of channels, 
    %                 T is Tf x nC x nS    
    % and the resulting 
    %                 M is (nS*nC) x Tf
    [Tf nC nS] = size(T);
    M = reshape(T, Tf, nS*nC)';