function VC = t2vc(T)
    % converts the tensor represented waveforms into vector concatenated
    % representation. T is a tensor of size Tf x nC x nS where
    % Tf is the number of samples per channel and waveform, nC is the
    % number of recording channels and nS is the number of waveforms
    
    [Tf, nC, nS] = size(T);    
    VC = reshape(permute(T, [2 1 3]), nC, Tf*nS);