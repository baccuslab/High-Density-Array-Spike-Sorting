function V = t2v(T)
    % converts the tensor representation of waveforms into the vector
    % representation. T is a tensort of size Tf x nC x nS where Tf is the
    % number of samples per waveform on one cahnnel, nC is the number of
    % channels and nS is the number of waveforms
    %
    % This also works if the subsample shifted versions of T are stored in
    % the 4th dimension. Then, V will be a tensor (nS x Tf*nC x nSumsampleTaus)
    [Tf, nC, nS, nSubsampleTaus] = size(T);
    P = reshape(T, Tf*nC, nS, nSubsampleTaus);
    V = permute(P, [2 1 3]);
